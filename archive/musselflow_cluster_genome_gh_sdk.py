"""
Decode a compact Galapagos genome into clustered star-control values.

Use one Gene Pool containing RequiredGenes values in the 0.0-1.0 domain.
Eight genes describe each cluster:
    U, V, spread, aspect, heading, rod-A length, rod-B length, twist.

    Inputs:
        run: Generate the field {item,bool}
        domain: Surface on which the stars are distributed {item,surface}
        count: Total number of stars {item,int}
        clusterCount: Number of interacting spatial families {item,int}
        genes: Normalized Gene Pool values, each 0.0-1.0 {list,float}
        flowDir: Reference water-flow direction {item,vector}
        minLength: Minimum rod length {item,float}
        maxLength: Maximum rod length {item,float}

    Outputs (this exact order):
        RodValues: Two dimensions per star, flattened as
                   [A0, B0, A1, B1, ...]. Partition with size 2.
        Centers: One center per star
        StarID: One star index per pair
        ClusterID: Cluster index per star
        LengthA: First rod length per star
        LengthB: Second rod length per star
        Angles: First rod angle in degrees relative to flow
        RequiredGenes: Required Gene Pool size
        Report: Genome and geometry summary
"""

import System
import Rhino
import Grasshopper
import math

rc = Rhino


def clamp(value, low=0.0, high=1.0):
    return max(low, min(high, float(value)))


def remap01(value, low, high):
    return low + clamp(value)*(high-low)


def normalized_genes(values):
    if values is None:
        return []
    try:
        return [clamp(value) for value in values]
    except TypeError:
        return [clamp(values)]


def get_surface(value):
    """Accept Surface, BrepFace, or a single-face Brep."""
    if value is None:
        return None
    if hasattr(value, "Value"):
        value = value.Value
    if isinstance(value, rc.Geometry.Surface):
        return value
    if isinstance(value, rc.Geometry.BrepFace):
        return value
    if isinstance(value, rc.Geometry.Brep) and value.Faces.Count:
        return value.Faces[0]
    return None


def surface_frame(surface, u01, v01):
    u = surface.Domain(0).ParameterAt(clamp(u01))
    v = surface.Domain(1).ParameterAt(clamp(v01))
    success, frame = surface.FrameAt(u, v)
    if success:
        return frame
    point = surface.PointAt(u, v)
    return rc.Geometry.Plane(point, rc.Geometry.Vector3d.ZAxis)


def vector_in_plane(vector, plane):
    """Project flow into the local surface plane."""
    result = rc.Geometry.Vector3d(vector)
    result -= plane.Normal * (result * plane.Normal)
    if result.Length < 1e-9:
        result = rc.Geometry.Vector3d(plane.XAxis)
    result.Unitize()
    return result


def rotated(vector, radians, axis):
    result = rc.Geometry.Vector3d(vector)
    result.Rotate(radians, axis)
    result.Unitize()
    return result


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def RunScript(
            self,
            run: bool,
            domain: Rhino.Geometry.Surface,
            count: int,
            clusterCount: int,
            genes: System.Collections.Generic.List[float],
            flowDir: Rhino.Geometry.Vector3d,
            minLength: float,
            maxLength: float):

        empty = ([], [], [], [], [], [], [], 0, ["Waiting for inputs."])
        if not run:
            return empty

        surface = get_surface(domain)
        if surface is None:
            return ([], [], [], [], [], [], [], 0,
                    ["Connect one valid Surface to domain."])

        star_count = max(1, int(count) if count else 20)
        cluster_count = max(1, min(
            star_count, int(clusterCount) if clusterCount else 3))
        required = cluster_count*8
        genome = normalized_genes(genes)
        if len(genome) < required:
            return ([], [], [], [], [], [], [], required, [
                "Gene Pool too small: received %d, requires %d."
                % (len(genome), required),
                "Set every gene domain to 0.0-1.0."])

        low_length = max(0.0, float(minLength) if minLength else 1.0)
        high_length = max(
            low_length, float(maxLength) if maxLength else low_length*4.0)
        flow = (rc.Geometry.Vector3d(flowDir)
                if flowDir is not None and flowDir.Length > 1e-9
                else rc.Geometry.Vector3d.XAxis)

        rod_values = []
        centers = []
        star_ids = []
        cluster_ids = []
        lengths_a = []
        lengths_b = []
        angles = []

        # Assign almost equal numbers of stars to each cluster.
        local_indices = [0]*cluster_count
        cluster_sizes = [
            star_count//cluster_count +
            (1 if cluster < star_count % cluster_count else 0)
            for cluster in range(cluster_count)]
        golden_angle = math.pi*(3.0-math.sqrt(5.0))

        for star in range(star_count):
            cluster = star % cluster_count
            local = local_indices[cluster]
            local_indices[cluster] += 1
            size = cluster_sizes[cluster]
            gene = genome[cluster*8:(cluster+1)*8]

            center_u = gene[0]
            center_v = gene[1]
            spread = remap01(gene[2], 0.025, 0.32)
            aspect = remap01(gene[3], 0.35, 1.65)
            heading = remap01(gene[4], -math.pi, math.pi)
            base_a = remap01(gene[5], low_length, high_length)
            base_b = remap01(gene[6], low_length, high_length)
            twist = remap01(gene[7], -math.pi, math.pi)

            # Deterministic sunflower distribution: stable fitness for Galapagos.
            radius = math.sqrt((local+0.5)/float(max(1, size)))
            phase = local*golden_angle
            du = spread*radius*math.cos(phase)*aspect
            dv = spread*radius*math.sin(phase)/aspect
            cos_h = math.cos(heading)
            sin_h = math.sin(heading)
            u01 = clamp(center_u + du*cos_h-dv*sin_h)
            v01 = clamp(center_v + du*sin_h+dv*cos_h)

            frame = surface_frame(surface, u01, v01)
            flow_axis = vector_in_plane(flow, frame)

            # Radial gradient makes members related but not identical.
            local_angle = heading + twist*radius
            modulation = 0.72 + 0.28*math.cos(phase-heading)
            length_a = max(low_length, base_a*modulation)
            length_b = max(
                low_length, base_b*(1.72-modulation))

            rod_values.extend([length_a, length_b])
            centers.append(frame.Origin)
            star_ids.append(star)
            cluster_ids.append(cluster)
            lengths_a.append(length_a)
            lengths_b.append(length_b)
            angles.append(math.degrees(local_angle))

        report = [
            "CLUSTER GENOME | %d stars | %d clusters | %d/%d genes"
            % (star_count, cluster_count, required, len(genome)),
            "RodValues contains %d values: consecutive pairs A,B for %d stars."
            % (len(rod_values), star_count),
            "Use Partition List with size 2; multiply each A and B four times "
            "in your existing eight-edge star construction.",
            "Genes are deterministic: the same genome always gives the same field.",
            "Keep Galapagos mutation low but nonzero; maximize the simulation "
            "Fitness output, not an animation frame."
        ]

        return (
            rod_values,
            centers,
            star_ids,
            cluster_ids,
            lengths_a,
            lengths_b,
            angles,
            required,
            report)
