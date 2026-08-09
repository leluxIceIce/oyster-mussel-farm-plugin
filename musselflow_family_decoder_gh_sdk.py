"""
Expand three small Gene Pools into the values required by an existing
Grasshopper star-field definition.

This component does not replace the downstream geometry algorithm. Insert it
between the existing red/yellow Gene Pools and the existing Evaluate Surface /
rod construction:

    U Gene Pool ----> uGenes      UValues ----> existing surface U input
    V Gene Pool ----> vGenes      VValues ----> existing surface V input
    Rod Gene Pool --> rodGenes    RodValues --> Partition List (size 2)

Recommended: 8 genes in each pool, all with domain 0.0 to 1.0.

    Inputs:
        run: Expand the control genes {item,bool}
        uGenes: Eight normalized U-family controls {list,float}
        vGenes: Eight normalized V-family controls {list,float}
        rodGenes: Eight normalized rod-family controls {list,float}
        count: Required number of stars, e.g. 20 {item,int}
        rodDomain: Physical output interval for rod dimensions {item,interval}
        restart: Button; advances the internal starting seed once {item,bool}

    Outputs (this exact order):
        UValues: One normalized U coordinate per star
        VValues: One normalized V coordinate per star
        RodValues: Flat [A0,B0,A1,B1,...] list; partition with size 2
        Angles: Optional orientation angle per star, in degrees
        FamilyID: Family/cluster index per star
        Report: Generator settings and validation
"""

import System
import Rhino
import Grasshopper
import math

rc = Rhino


def clamp(value, low=0.0, high=1.0):
    return max(low, min(high, float(value)))


def pool(values):
    if values is None:
        return []
    try:
        return [clamp(value) for value in values]
    except TypeError:
        return [clamp(values)]


def remap(value, low, high):
    return low + clamp(value)*(high-low)


def fract(value):
    return value-math.floor(value)


def edge_spread(value, exponent):
    """Monotonic contrast curve that moves values toward 0 and 1."""
    value = clamp(value)
    if value <= 0.5:
        return 0.5*math.pow(2.0*value, exponent)
    return 1.0-0.5*math.pow(2.0*(1.0-value), exponent)


def normalize(values):
    low = min(values)
    high = max(values)
    span = high-low
    if span < 1e-12:
        return [0.5]*len(values)
    return [(value-low)/span for value in values]


def lorenz_sequence(count, ug, vg, rg, seed):
    """Bounded deterministic samples of a Lorenz trajectory."""
    sigma = remap(ug[5], 8.0, 14.0)
    rho = remap(vg[5], 22.0, 34.0)
    beta = remap(rg[7], 2.2, 3.2)
    seed_u = fract(seed*0.61803398875)
    seed_v = fract(seed*0.41421356237)
    seed_z = fract(seed*0.73205080757)
    x = remap(fract(ug[0]+seed_u), -1.0, 1.0)
    y = remap(fract(vg[0]+seed_v), -1.0, 1.0)
    z = remap(fract(rg[0]+seed_z), 0.5, 2.0)
    dt = 0.008
    burn = 80
    stride = 5
    xs, ys, zs = [], [], []
    total = burn + count*stride
    for step in range(total):
        dx = sigma*(y-x)
        dy = x*(rho-z)-y
        dz = x*y-beta*z
        x += dx*dt
        y += dy*dt
        z += dz*dt
        if step >= burn and (step-burn) % stride == 0:
            xs.append(x)
            ys.append(y)
            zs.append(z)
    return normalize(xs), normalize(ys), normalize(zs)


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def RunScript(
            self,
            run: bool,
            uGenes: System.Collections.Generic.List[float],
            vGenes: System.Collections.Generic.List[float],
            rodGenes: System.Collections.Generic.List[float],
            count: int,
            rodDomain: Rhino.Geometry.Interval,
            restart: bool):

        if not hasattr(self, "_internalSeed"):
            self._internalSeed = 0
            self._restartWasPressed = False

        pressed = bool(restart)
        if pressed and not self._restartWasPressed:
            self._internalSeed += 1
        self._restartWasPressed = pressed
        seed = self._internalSeed

        empty = ([], [], [], [], [], [
            "Waiting for inputs. Internal seed: %d." % seed])
        if not run:
            return empty

        ug = pool(uGenes)
        vg = pool(vGenes)
        rg = pool(rodGenes)
        if len(ug) < 8 or len(vg) < 8 or len(rg) < 8:
            return ([], [], [], [], [], [
                "Use at least 8 genes in each of the U, V and Rod Gene Pools.",
                "Set every Gene Pool domain to 0.0-1.0."])

        star_count = max(1, int(count) if count else 20)
        if rodDomain is not None and rodDomain.IsValid:
            rod_low = min(rodDomain.T0, rodDomain.T1)
            rod_high = max(rodDomain.T0, rodDomain.T1)
        else:
            rod_low, rod_high = 0.7, 1.2

        lx, ly, lz = lorenz_sequence(star_count, ug, vg, rg, seed)

        # Family count, centres and spread are controlled by the red pools.
        family_count = 2 + int(round(4.0*0.5*(ug[4]+vg[4])))
        global_u = remap(ug[0], 0.10, 0.90)
        global_v = remap(vg[0], 0.10, 0.90)
        center_radius = remap(0.5*(ug[2]+vg[2]), 0.05, 0.38)
        spread = remap(0.5*(ug[1]+vg[1]), 0.03, 0.24)
        aspect = remap(0.5*(ug[3]+vg[3]), 0.45, 2.10)
        phase = 2.0*math.pi*fract(
            0.5*(ug[7]+vg[7]) + seed*0.61803398875)
        uv_edge_exponent = remap(
            0.5*(ug[5]+vg[5]), 1.25, 1.85)

        # Smoothly blend clusters, a directional gradient and the attractor.
        attractor_blend = remap(0.5*(ug[6]+vg[6]), 0.0, 0.55)
        gradient_blend = remap(0.5*(ug[7]+vg[7]), 0.0, 0.30)
        cluster_blend = max(0.0, 1.0-attractor_blend-gradient_blend)
        golden_angle = math.pi*(3.0-math.sqrt(5.0))

        u_values = []
        v_values = []
        rod_raw_values = []
        rod_values = []
        angles = []
        family_ids = []

        for index in range(star_count):
            family = index % family_count
            local = index // family_count
            family_size = (
                star_count//family_count +
                (1 if family < star_count % family_count else 0))
            radial = math.sqrt((local+0.5)/float(max(1, family_size)))
            family_phase = phase + 2.0*math.pi*family/family_count
            local_phase = golden_angle*local + family_phase

            family_u = (
                global_u + center_radius*math.cos(family_phase)
                + spread*aspect*radial*math.cos(local_phase))
            family_v = (
                global_v + center_radius*math.sin(family_phase)
                + spread/aspect*radial*math.sin(local_phase))

            t = (index+0.5)/float(star_count)
            gradient_u = clamp(
                t + (ug[3]-0.5)*0.35*math.sin(2.0*math.pi*t+phase))
            gradient_v = clamp(
                global_v + (vg[3]-0.5)*(2.0*t-1.0)
                + (vg[2]-0.5)*0.25*math.sin(4.0*math.pi*t+phase))

            u = edge_spread(clamp(
                cluster_blend*family_u +
                gradient_blend*gradient_u +
                attractor_blend*lx[index]), uv_edge_exponent)
            v = edge_spread(clamp(
                cluster_blend*family_v +
                gradient_blend*gradient_v +
                attractor_blend*ly[index]), uv_edge_exponent)

            # Two related dimensions per star; downstream repeats each four times.
            family_wave = 0.5+0.5*math.sin(
                family_phase + radial*math.pi)
            contrast = remap(rg[6], 0.30, 1.00)
            length_a_01 = clamp(
                rg[0] + contrast*(
                    0.32*(2.0*family_wave-1.0)
                    + 0.22*(2.0*t-1.0)
                    + 0.18*(2.0*lz[index]-1.0)
                    + 0.12*(2.0*rg[2]-1.0)))
            length_b_01 = clamp(
                rg[1] + contrast*(
                    0.32*(1.0-2.0*family_wave)
                    + 0.22*(1.0-2.0*t)
                    + 0.18*(2.0*lz[index]-1.0)
                    + 0.12*(2.0*rg[3]-1.0)))

            u_values.append(u)
            v_values.append(v)
            rod_raw_values.extend([length_a_01, length_b_01])
            family_ids.append(family)

        # Expand the actual population range before mapping to physical lengths.
        # This retains each pattern's ordering but prevents middle-heavy rods.
        raw_low = min(rod_raw_values)
        raw_high = max(rod_raw_values)
        raw_span = raw_high-raw_low
        range_mix = remap(rg[7], 0.72, 0.94)
        rod_edge_exponent = remap(rg[6], 1.15, 1.65)
        for value in rod_raw_values:
            stretched = ((value-raw_low)/raw_span
                         if raw_span > 1e-12 else 0.5)
            expanded = (1.0-range_mix)*value + range_mix*stretched
            expanded = edge_spread(expanded, rod_edge_exponent)
            rod_values.append(remap(expanded, rod_low, rod_high))

        # Tangent of the generated UV trajectory gives an optional orientation.
        for index in range(star_count):
            previous = (index-1) % star_count
            following = (index+1) % star_count
            du = u_values[following]-u_values[previous]
            dv = v_values[following]-v_values[previous]
            base_angle = math.degrees(math.atan2(dv, du))
            alternating = remap(rg[6], -45.0, 45.0)
            angles.append(base_angle + alternating*(1 if index % 2 else -1))

        report = [
            "FAMILY DECODER | %d stars | %d families | 24 control genes | "
            "seed %d" % (star_count, family_count, seed),
            "Outputs: %d U + %d V + %d rod values (%d A/B pairs)."
            % (len(u_values), len(v_values), len(rod_values), star_count),
            "Blend: clusters %.2f | gradient %.2f | Lorenz-inspired %.2f"
            % (cluster_blend, gradient_blend, attractor_blend),
            "Ranges: U %.3f-%.3f | V %.3f-%.3f | rods %.3f-%.3f"
            % (min(u_values), max(u_values), min(v_values), max(v_values),
               min(rod_values), max(rod_values)),
            "The attractor generates geometry only; it is not the water solver.",
            "Same genes always produce the same values, so Galapagos fitness "
            "remains deterministic until restart is pressed.",
            "Press restart only before starting a new Galapagos run."
        ]

        return (
            u_values,
            v_values,
            rod_values,
            angles,
            family_ids,
            report)
