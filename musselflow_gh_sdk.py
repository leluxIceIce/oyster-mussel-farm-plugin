# r: numpy
"""
Plan-view flow analysis of water passing a formation of mussel dropper-lines.
SDK-mode version: typed RunScript inputs/outputs (the node self-declares its
type hints) plus a custom viewport preview that draws a color legend, a flow
arrow and the domain frame -- so the analysis reads like a chart inside Rhino,
no baking or PNG needed.

Looks straight down at the water: water flows in FlowDir, each dropper line
(buoy line + mesh sock + mussels) is a vertical body seen as a cross-section from
above. Socks may be Points (circle via Radius), closed Curves, Meshes or
Surfaces; 3D socks are projected into plan silhouettes so every supplied object
participates even when geometries occupy different Z levels. Numbers are
dimensionless -- use them to RANK layouts.

Paste into a Python 3 Script component, then switch it to SDK mode (dashboard:
"Convert to SDK mode" / the class below). numpy loads via the `# r: numpy` line.
-
Name: MusselFlowAnalysisSDK
Updated: 260801
Author: Felix (with Claude)
Copyright: Creative Commons - Attribution 4.0 International

    Inputs (RunScript annotations set these hints automatically):
        run: Activates the component {item,bool}
        socks: Sock geometry -- Points, closed Curves, Meshes or Surfaces;
               every Grasshopper tree branch is flattened into one farm
               {tree,geometry}
        radius: Sock radius in model units, for Point socks {item,float}
        flowDir: Concurrent ambient-flow samples; all branches are flattened
                 and averaged into one representative vector {tree,vector}
        domain: Optional Rectangle bounding the area, else auto {item,rectangle}
        resolution: Grid cells along the flow axis, 120-260 {item,int}
        regions: Optional closed Curves = measurement patches {list,curve}
        uptake: Mussel filtering strength 0.0-0.5 {item,float}
    Outputs:
        FlowMesh: Mesh colored by flow slowdown for the selected scenario
        FoodMesh: Mesh colored by seston/food for the selected scenario
        Report: One text line per region with its scores
        Delivery: Per-region mean flow delivery (1.0 = full freestream)
        Deflection: Per-region mean flow deflection, degrees
        Filtered: Per-region mean fraction of seston removed

The legend, colour scales, flow arrow and panel frames are drawn automatically in
the viewport (DrawViewportWires/Meshes) whenever the component previews -- you do
NOT need to wire any output to see them.
"""

import System
import Rhino
import Grasshopper
import numpy as np
import math
import time

rc = Rhino


# ================================================================ numpy core
E = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
              [1, 1], [-1, 1], [-1, -1], [1, -1]], dtype=int)
W = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])


def equilibrium(rho, ux, uy):
    """ D2Q9 lattice-Boltzmann equilibrium distribution """
    usq = ux*ux + uy*uy
    feq = np.empty((9,) + rho.shape)
    for i in range(9):
        cu = E[i, 0]*ux + E[i, 1]*uy
        feq[i] = W[i]*rho*(1 + 3*cu + 4.5*cu*cu - 1.5*usq)
    return feq


def solveFlow(nx, ny, solid, u0, tau, steps, tol=1e-6, check=100):
    """ Lattice-Boltzmann velocity field. Flow enters +x left, leaves right,
    periodic in y (one row = an infinite line of droppers). Returns ux, uy, speed. """
    omega = 1.0/tau
    ux = np.full((ny, nx), u0)
    uy = np.zeros((ny, nx))
    rho = np.ones((ny, nx))
    f = equilibrium(rho, ux, uy)
    fluid = ~solid
    prev = None
    for t in range(steps):
        f[[3, 6, 7], :, -1] = f[[3, 6, 7], :, -2]
        rho = f.sum(0)
        ux = (f*E[:, 0][:, None, None]).sum(0)/rho
        uy = (f*E[:, 1][:, None, None]).sum(0)/rho
        ux[:, 0] = u0
        uy[:, 0] = 0.0
        rho[:, 0] = 1.0
        feq = equilibrium(rho, ux, uy)
        f[:, :, 0] = feq[:, :, 0]
        fout = f - omega*(f - feq)
        for i in range(9):
            fout[i, solid] = f[OPP[i], solid]
        for i in range(9):
            f[i] = np.roll(fout[i], (E[i, 1], E[i, 0]), axis=(0, 1))
        if t % check == 0:
            sp = np.sqrt(ux*ux + uy*uy)*fluid
            if prev is not None and np.abs(sp - prev).mean() < tol:
                break
            prev = sp
    ux *= fluid
    uy *= fluid
    return ux, uy, np.sqrt(ux*ux + uy*uy)


def solveSeston(nx, ny, ux, uy, solid, band, uptake, D, steps=4000, tol=1e-7):
    """ Advect-diffuse seston with an uptake sink in the mussel band. FTCS, so
    D <= 0.25 for stability. Returns the seston field. """
    D = min(D, 0.24)
    s = np.ones((ny, nx))
    fluid = ~solid
    for t in range(steps):
        s[:, 0] = 1.0
        adv = (np.maximum(ux, 0)*(s - np.roll(s, 1, 1)) +
               np.minimum(ux, 0)*(np.roll(s, -1, 1) - s) +
               np.maximum(uy, 0)*(s - np.roll(s, 1, 0)) +
               np.minimum(uy, 0)*(np.roll(s, -1, 0) - s))
        lap = (np.roll(s, 1, 1) + np.roll(s, -1, 1) +
               np.roll(s, 1, 0) + np.roll(s, -1, 0) - 4*s)
        sn = s - adv + D*lap
        sn[band] *= (1 - uptake)
        sn[solid] = 0.0
        sn = np.clip(sn, 0.0, 1.5)
        sn[:, 0] = 1.0
        if t % 100 == 0 and t > 0 and np.abs(sn - s).mean() < tol:
            s = sn
            break
        s = sn
    return s*fluid


def rasterizePolygons(nx, ny, polys, ensure_each=False):
    """Fill cells whose *centres* fall inside any footprint polygon.

    The previous version tested integer grid corners although ``cellToWorld``
    places mesh vertices at ``i + 0.5, j + 0.5``. Narrow socks could therefore
    be accepted as geometry yet vanish from the solid mask. Work is now limited
    to each polygon's bounding box, and ``ensure_each`` conservatively assigns
    a valid sub-cell footprint to its nearest cell.
    """
    inside = np.zeros((ny, nx), dtype=bool)
    for raw_poly in polys:
        poly = np.asarray(raw_poly, dtype=float)
        if len(poly) < 3 or not np.isfinite(poly).all():
            continue

        min_x = max(0, int(math.floor(float(poly[:, 0].min())-0.5)))
        max_x = min(nx-1, int(math.ceil(float(poly[:, 0].max())-0.5)))
        min_y = max(0, int(math.floor(float(poly[:, 1].min())-0.5)))
        max_y = min(ny-1, int(math.ceil(float(poly[:, 1].max())-0.5)))
        polygon_mask = np.zeros((ny, nx), dtype=bool)

        if min_x <= max_x and min_y <= max_y:
            ys, xs = np.mgrid[min_y:max_y+1, min_x:max_x+1]
            px = xs.ravel().astype(float)+0.5
            py = ys.ravel().astype(float)+0.5
            vx, vy = poly[:, 0], poly[:, 1]
            local = np.zeros(px.shape, dtype=bool)
            previous = len(poly)-1
            for current in range(len(poly)):
                crosses = ((vy[current] > py) != (vy[previous] > py)) & (
                    px < (vx[previous]-vx[current])*(py-vy[current]) /
                    (vy[previous]-vy[current]+1e-12)+vx[current])
                local ^= crosses
                previous = current
            polygon_mask[min_y:max_y+1, min_x:max_x+1] = local.reshape(
                max_y-min_y+1, max_x-min_x+1)

        if ensure_each and not polygon_mask.any():
            # A real obstruction smaller than one cell must still participate.
            center = poly.mean(axis=0)
            cell_x = int(round(float(center[0])-0.5))
            cell_y = int(round(float(center[1])-0.5))
            if 0 <= cell_x < nx and 0 <= cell_y < ny:
                polygon_mask[cell_y, cell_x] = True
        inside |= polygon_mask
    return inside


def bandAround(mask, width):
    """ Mussel filter band: fluid cells within 'width' of the sock surface, any
    shape. Dilation(mask) minus mask. """
    out = mask.copy()
    for _ in range(int(round(width))):
        out = (out | np.roll(out, 1, 0) | np.roll(out, -1, 0) |
               np.roll(out, 1, 1) | np.roll(out, -1, 1))
    return out & ~mask


# ============================================================= color mapping
def lerpColor(stops, t):
    """ Interpolate through (position,(r,g,b)) color stops """
    t = 0.0 if t < 0 else 1.0 if t > 1 else t
    for k in range(len(stops) - 1):
        t0, c0 = stops[k]
        t1, c1 = stops[k + 1]
        if t <= t1:
            f = 0.0 if t1 == t0 else (t - t0)/(t1 - t0)
            return tuple(int(c0[i] + f*(c1[i] - c0[i])) for i in range(3))
    return stops[-1][1]


def cmapDeflect(t):
    """ blue undisturbed -> cyan -> yellow -> orange -> red -> purple stopped """
    stops = [(0.0, (0, 51, 204)), (0.25, (0, 179, 179)), (0.5, (255, 224, 0)),
             (0.7, (255, 136, 0)), (0.88, (224, 0, 0)), (1.0, (122, 0, 160))]
    return lerpColor(stops, t)


def cmapFood(t):
    """ dark starved -> teal -> yellow full food """
    stops = [(0.0, (40, 0, 60)), (0.5, (33, 145, 140)), (1.0, (253, 231, 37))]
    return lerpColor(stops, t)


# ================================================= Rhino geometry -> polygons
def toList(x):
    """ Coerce any GH input into a python list: passes lists through, wraps a
    single item, and treats a lone Guid/None safely (a Guid is not iterable, which
    is what crashes list(socks) when the input is Item access or ghdoc-hinted). """
    if x is None:
        return []
    if isinstance(x, System.Guid):
        return [x]
    if hasattr(x, "Branches"):
        try:
            return [item for branch in x.Branches for item in branch]
        except Exception:
            pass
    try:
        return list(x)
    except TypeError:
        return [x]


def flowVectorList(value):
    """Resolve one/list GH vector input into finite nonzero Vector3d values."""
    vectors = []
    for candidate in toList(value):
        if hasattr(candidate, "Value"):
            candidate = candidate.Value
        if not isinstance(candidate, rc.Geometry.Vector3d):
            continue
        try:
            valid = candidate.IsValid and candidate.Length > 1e-9
            finite = all(math.isfinite(component) for component in (
                candidate.X, candidate.Y, candidate.Z))
        except Exception:
            valid = finite = False
        if valid and finite:
            vectors.append(candidate)
    return vectors


def resolveGeo(o):
    """ Return real RhinoCommon geometry from a GH input item, whatever form it
    arrives in: unwrap GH_ wrappers (.Value), and resolve document GUIDs (the
    'ghdoc Object' type hint marshals geometry to Guids) back to their geometry. """
    # GH_ObjectWrapper and typed GH_Goo can be nested. Unwrap repeatedly, but
    # stop when a wrapper returns itself.
    for _ in range(4):
        if not hasattr(o, "Value"):
            break
        value = o.Value
        if value is o:
            break
        o = value
    if hasattr(o, "ScriptVariable"):
        try:
            value = o.ScriptVariable()
            if value is not None:
                o = value
        except Exception:
            pass
    if isinstance(o, System.Guid):
        obj = None
        try:
            import scriptcontext as sc
            obj = sc.doc.Objects.FindId(o)
        except Exception:
            obj = None
        if obj is None:
            doc = rc.RhinoDoc.ActiveDoc
            if doc is not None:
                obj = doc.Objects.FindId(o)
        if obj is not None:
            return obj.Geometry
    return o


def unwrap(o):
    """ Back-compat shim: resolve a single item to geometry """
    return resolveGeo(o)


def meanGeoZ(geos):
    """ Mean bounding-box-center Z of the sock geometry (fallback analysis height) """
    zs = []
    for g in geos:
        if isinstance(g, rc.Geometry.Point3d):
            zs.append(g.Z)
        elif isinstance(g, rc.Geometry.Point):
            zs.append(g.Location.Z)
        else:
            try:
                zs.append(g.GetBoundingBox(True).Center.Z)
            except Exception:
                pass
    return sum(zs)/len(zs) if zs else 0.0


def domainInfo(domain):
    """ Resolve a domain input into (analysisZ or None, list of world (x,y) boundary
    points or None). Accepts a Rectangle3d, a rectangle Curve, a document Guid (the
    ghdoc hint turns a rectangle into a curve), or None. """
    if domain is None:
        return None, None
    d = resolveGeo(domain)
    G = rc.Geometry
    if isinstance(d, G.Rectangle3d):
        pts = [(d.Corner(i).X, d.Corner(i).Y) for i in range(4)]
        return d.Plane.OriginZ, pts
    if isinstance(d, G.Curve):
        try:
            z = d.GetBoundingBox(True).Center.Z
        except Exception:
            z = None
        return z, curvePoints(d, 64)
    return None, None


def curvePoints(crv, n=90):
    """ Sample a curve into a list of (x,y) tuples (its footprint outline) """
    ts = crv.DivideByCount(n, True)
    if not ts:
        p0, p1 = crv.PointAtStart, crv.PointAtEnd
        return [(p0.X, p0.Y), (p1.X, p1.Y)]
    return [(crv.PointAt(t).X, crv.PointAt(t).Y) for t in ts]


def circlePoints(cx, cy, r, k=48):
    """ A circle outline in world XY as (x,y) tuples """
    return [(cx + r*math.cos(a), cy + r*math.sin(a))
            for a in np.linspace(0, 2*math.pi, k, endpoint=False)]


def convexHull2D(points):
    """Return a deterministic fallback hull for projected mesh vertices."""
    unique = sorted(set((float(point[0]), float(point[1])) for point in points))
    if len(unique) <= 2:
        return unique

    def cross(origin, first, second):
        return ((first[0]-origin[0])*(second[1]-origin[1]) -
                (first[1]-origin[1])*(second[0]-origin[0]))

    lower = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1]+upper[:-1]


def projectedMeshPolys(mesh, plane):
    """Obtain plan silhouettes, with a vertex hull as a robust fallback."""
    polygons = []
    try:
        outlines = mesh.GetOutlines(plane)
    except Exception:
        outlines = None
    if outlines:
        for outline in outlines:
            points = [(point.X, point.Y) for point in outline]
            if len(points) >= 3:
                polygons.append(points)
    if polygons:
        return polygons
    vertices = [
        (mesh.Vertices[index].X, mesh.Vertices[index].Y)
        for index in range(mesh.Vertices.Count)]
    hull = convexHull2D(vertices)
    return [hull] if len(hull) >= 3 else []


def brepProjectedPolys(brep, plane):
    """Mesh a Brep for plan-outline extraction without changing the source."""
    try:
        meshes = rc.Geometry.Mesh.CreateFromBrep(
            brep, rc.Geometry.MeshingParameters.FastRenderMesh)
    except Exception:
        meshes = None
    if not meshes:
        points = [(vertex.Location.X, vertex.Location.Y)
                  for vertex in brep.Vertices]
        hull = convexHull2D(points)
        return [hull] if len(hull) >= 3 else []
    joined = rc.Geometry.Mesh()
    for mesh in meshes:
        if mesh is not None:
            joined.Append(mesh)
    return projectedMeshPolys(joined, plane)


def sockToWorldPolys(geo, radius, plane, tol):
    """Convert each sock into a plan-projected world-space obstruction."""
    G = rc.Geometry
    geo = unwrap(geo)
    out = []
    if isinstance(geo, G.Point3d):
        out.append(circlePoints(geo.X, geo.Y, radius))
    elif isinstance(geo, G.Point):
        out.append(circlePoints(geo.Location.X, geo.Location.Y, radius))
    elif isinstance(geo, G.Curve):
        out.append(curvePoints(geo))
    elif isinstance(geo, G.Mesh):
        out.extend(projectedMeshPolys(geo, plane))
    elif isinstance(geo, G.Extrusion):
        out.extend(brepProjectedPolys(geo.ToBrep(), plane))
    elif isinstance(geo, G.Brep):
        out.extend(brepProjectedPolys(geo, plane))
    elif isinstance(geo, G.Surface):
        out.extend(brepProjectedPolys(geo.ToBrep(), plane))
    elif hasattr(G, "SubD") and isinstance(geo, G.SubD):
        try:
            out.extend(brepProjectedPolys(geo.ToBrep(), plane))
        except Exception:
            pass
    return out


def buildGrid(boundPts, flowDir, res, radius, domainPts):
    """ Flow-aligned solve grid (xhat along flow, yhat 90deg left). Returns sizing
    + world<->grid mappers + frame vectors. boundPts: list of (x,y) to enclose.
    domainPts: optional list of world (x,y) that override the auto extents. """
    fx, fy = flowDir
    nrm = math.hypot(fx, fy)
    xhat = (1.0, 0.0) if nrm < 1e-9 else (fx/nrm, fy/nrm)
    yhat = (-xhat[1], xhat[0])

    def proj(p, ax):
        return p[0]*ax[0] + p[1]*ax[1]

    us = [proj(p, xhat) for p in boundPts]
    vs = [proj(p, yhat) for p in boundPts]
    umin, umax, vmin, vmax = min(us), max(us), min(vs), max(vs)

    usePts = None                          # accept only a real list of (x,y) points;
    if domainPts and not isinstance(domainPts, System.Guid):
        try:                               # ignore a stray Guid/object so a domain
            usePts = [(float(p[0]), float(p[1])) for p in domainPts]
        except (TypeError, IndexError):    # never crashes the solve -> auto-domain
            usePts = None
    if usePts:
        cu = [proj(p, xhat) for p in usePts]
        cv = [proj(p, yhat) for p in usePts]
        umin, umax, vmin, vmax = min(cu), max(cu), min(cv), max(cv)
    else:
        spanU = max(umax - umin, 2*radius)
        spanV = max(vmax - vmin, 2*radius)
        umin -= 0.5*spanV + 1.5*spanU
        umax += 1.5*spanV + 4.0*spanU
        vmin -= 0.3*spanV + 1.5*spanU
        vmax += 0.3*spanV + 1.5*spanU

    Wm, Hm = umax - umin, vmax - vmin
    nx = int(res)
    h = Wm/nx
    ny = max(8, int(round(Hm/h)))
    ox = xhat[0]*umin + yhat[0]*vmin
    oy = xhat[1]*umin + yhat[1]*vmin

    def worldToGrid(px, py):
        dx, dy = px - ox, py - oy
        return (dx*xhat[0] + dy*xhat[1])/h, (dx*yhat[0] + dy*yhat[1])/h

    def cellToWorld(i, j):
        wx = ox + xhat[0]*((i + 0.5)*h) + yhat[0]*((j + 0.5)*h)
        wy = oy + xhat[1]*((i + 0.5)*h) + yhat[1]*((j + 0.5)*h)
        return wx, wy

    return dict(nx=nx, ny=ny, h=h, Wm=Wm, Hm=Hm, xhat=xhat, yhat=yhat,
                ox=ox, oy=oy,
                worldToGrid=worldToGrid, cellToWorld=cellToWorld)


def worldPolysToGrid(worldPolys, grid):
    """ Map world (x,y) outlines into (N,2) numpy arrays in grid coordinates """
    return [np.array([grid["worldToGrid"](x, y) for (x, y) in poly], float)
            for poly in worldPolys]


def domainMask(grid, domainPoints):
    """Return cells inside the supplied domain, or every cell for auto-domain."""
    if not domainPoints:
        return np.ones((grid["ny"], grid["nx"]), dtype=bool)
    polygons = worldPolysToGrid([domainPoints], grid)
    return rasterizePolygons(grid["nx"], grid["ny"], polygons)


def domainFrame(domainPoints, z, offset=(0.0, 0.0)):
    """Create the exact supplied-domain preview frame, optionally translated."""
    if not domainPoints:
        return None
    points = [
        rc.Geometry.Point3d(
            point[0]+offset[0], point[1]+offset[1], z)
        for point in domainPoints]
    if points and points[0].DistanceTo(points[-1]) > 1e-9:
        points.append(points[0])
    return rc.Geometry.Polyline(points)


def coloredMesh(field, cmap, solid, grid, z):
    """ Vertex-colored quad mesh over fluid cells, colored by 'field' via 'cmap' """
    mesh = rc.Geometry.Mesh()
    ny, nx = field.shape
    idx = np.full((ny, nx), -1, int)
    k = 0
    for j in range(ny):
        for i in range(nx):
            if solid[j, i]:
                continue
            wx, wy = grid["cellToWorld"](i, j)
            mesh.Vertices.Add(wx, wy, z)
            r, g, b = cmap(float(field[j, i]))
            mesh.VertexColors.Add(r, g, b)
            idx[j, i] = k
            k += 1
    for j in range(ny - 1):
        for i in range(nx - 1):
            a, b2, c, d = idx[j, i], idx[j, i+1], idx[j+1, i+1], idx[j+1, i]
            if a >= 0 and b2 >= 0 and c >= 0 and d >= 0:
                mesh.Faces.AddFace(int(a), int(b2), int(c), int(d))
    return mesh


def legendBar(ox, oy, z, xh, yh, width, height, cmap, segs=40):
    """ Gradient bar mesh rising along yh (value 0 -> 1), width along xh. Returns
    the mesh; colored by cmap so DrawMeshFalseColors shows the scale. """
    m = rc.Geometry.Mesh()
    for kk in range(segs + 1):
        t = kk/float(segs)
        px = ox + yh[0]*(t*height)
        py = oy + yh[1]*(t*height)
        m.Vertices.Add(px, py, z)
        m.Vertices.Add(px + xh[0]*width, py + xh[1]*width, z)
        r, g, b = cmap(t)
        m.VertexColors.Add(r, g, b)
        m.VertexColors.Add(r, g, b)
    for kk in range(segs):
        a = 2*kk
        m.Faces.AddFace(a, a+1, a+3, a+2)
    return m


def panelOverlay(corners, cmap, title, valueLabels, z, Wm, Hm, xh, yh, csize):
    """ Build the drawn furniture for one field panel: its frame polyline, its
    gradient legend bar (just downstream of the panel), and the text labels
    (title + value ticks). corners = (c00,c10,c11,c01) world (x,y) tuples with
    c00=upstream-bottom, c10=downstream-bottom, c11=downstream-top, c01=upstream-top. """
    c00, c10, c11, c01 = corners

    def P(xy, dx=0.0, dy=0.0):
        return rc.Geometry.Point3d(xy[0]+dx, xy[1]+dy, z)

    frame = rc.Geometry.Polyline([P(c00), P(c10), P(c11), P(c01), P(c00)])

    # legend bar sits just past the downstream (+xh) edge, rising along +yh
    bx = c10[0] + xh[0]*(0.05*Wm)
    by = c10[1] + xh[1]*(0.05*Wm)
    legW, legH = 0.04*Wm, 0.5*Hm
    bar = legendBar(bx, by, z, xh, yh, legW, legH, cmap)

    texts = [(title, P(c01, yh[0]*0.03*Hm, yh[1]*0.03*Hm))]
    for txt, val in valueLabels:
        lx = bx + yh[0]*(val*legH) + xh[0]*(legW + 0.3*csize)
        ly = by + yh[1]*(val*legH) + xh[1]*(legW + 0.3*csize)
        texts.append((txt, rc.Geometry.Point3d(lx, ly, z)))
    return bar, frame, texts


# ==================================================== SDK-mode GH component
class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    # ----- fixed solver settings (hidden to keep the node clean) -----
    U0 = 0.08
    TAU = 0.6
    # Upper bound only. The actual preview budget scales with grid length below;
    # a static 6000-step budget made selected-design inspection unnecessarily
    # slow and did not materially change the settled field in regression checks.
    STEPS = 1800
    BAND_W = 2
    SESTON_D = 0.15

    def RunScript(self,
                  run: bool,
                  socks: Grasshopper.DataTree[object],
                  radius: float,
                  flowDir: Grasshopper.DataTree[Rhino.Geometry.Vector3d],
                  domain: "Rhino.Geometry.Rectangle3d",
                  resolution: int,
                  regions: "System.Collections.Generic.List[Rhino.Geometry.Curve]",
                  uptake: float):
        """Average all valid vectors, then run the established solver once."""
        self._draw = None
        sockItems = toList(socks)
        if not run:
            return (None, None, [], [], [], [])
        if len(sockItems) == 0:
            return (None, None, [], [], [], [])

        vectors = flowVectorList(flowDir)
        used_default = not vectors
        if used_default:
            vectors = [rc.Geometry.Vector3d.XAxis]
        compositeWarning = None
        vector = rc.Geometry.Vector3d(
            sum(item.X for item in vectors)/len(vectors),
            sum(item.Y for item in vectors)/len(vectors),
            sum(item.Z for item in vectors)/len(vectors))
        if not vector.IsValid or vector.Length <= 1e-9:
            vector = vectors[0]
            compositeWarning = (
                "FLOW WARNING | Vector mean was near zero; first valid "
                "vector was used. Check for cancelling directions.")
        meanSpeed = sum(item.Length for item in vectors)/len(vectors)
        speedSpread = math.sqrt(sum(
            (item.Length-meanSpeed)**2 for item in vectors)/len(vectors))
        flowLabel = (
            "FLOW UNIFIED | %d vectors | mean (%.6f, %.6f, %.6f) | "
            "resultant %.6f | input-speed SD %.6f"
            % (len(vectors), vector.X, vector.Y, vector.Z,
               vector.Length, speedSpread))
        result = self._runSingle(
            run, sockItems, radius, vector, domain,
            resolution, regions, uptake)
        report = list(result[2])
        report.insert(0, flowLabel)
        report.insert(
            1, "FLOW MODEL LIMIT | Vector length is reported in input units, "
            "but this qualitative LBM uses fixed lattice inflow U0 = %.6f."
            % self.U0)
        if compositeWarning:
            report.insert(0, compositeWarning)
        if used_default:
            report.insert(
                0, "FLOW WARNING | No valid flow vectors; World X was used.")
        return (result[0], result[1], report,
                result[3], result[4], result[5])

    def _runSingle(self,
                   run: bool,
                   socks,
                   radius: float,
                   flowDir: "Rhino.Geometry.Vector3d",
                   domain: "Rhino.Geometry.Rectangle3d",
                   resolution: int,
                   regions: "System.Collections.Generic.List[Rhino.Geometry.Curve]",
                   uptake: float):
        """Run the established solver once for one unified flow vector."""
        tAll = time.perf_counter()
        self._draw = None
        FlowMesh = FoodMesh = None
        Report, Delivery, Deflection, Filtered = [], [], [], []

        if not run or len(toList(socks)) == 0:
            return (FlowMesh, FoodMesh, Report,
                    Delivery, Deflection, Filtered)

        suppliedItems = toList(socks)
        resolvedPairs = []
        unresolvedIndices = []
        for index, item in enumerate(suppliedItems):
            geometry = resolveGeo(item)
            if geometry is None:
                unresolvedIndices.append(index)
            else:
                resolvedPairs.append((index, geometry))
        sockGeos = [geometry for _, geometry in resolvedPairs]
        regionCrvs = [resolveGeo(c) for c in toList(regions)]
        regionCrvs = [c for c in regionCrvs
                      if isinstance(c, rc.Geometry.Curve)]
        if not sockGeos:
            return (FlowMesh, FoodMesh,
                    ["No sock geometry resolved -- check the Socks input."],
                    Delivery, Deflection, Filtered)

        tGeometry = time.perf_counter()
        doc = rc.RhinoDoc.ActiveDoc
        tol = doc.ModelAbsoluteTolerance if doc else 0.001
        domainZv, domainPts = domainInfo(domain)
        Z = domainZv if domainZv is not None else meanGeoZ(sockGeos)
        plane = rc.Geometry.Plane(rc.Geometry.Point3d(0, 0, Z),
                                  rc.Geometry.Vector3d.ZAxis)
        flow = ((flowDir.X, flowDir.Y)
                if (flowDir and flowDir.Length > 1e-9) else (1.0, 0.0))
        res = int(resolution) if resolution else 180
        res = max(80, min(res, 320))
        upt = uptake if uptake is not None else 0.15
        rad = radius if radius else 0.0

        # --- sock footprints (auto radius if missing) ---
        footprintGroups = []
        rejectedIndices = list(unresolvedIndices)
        for index, geometry in resolvedPairs:
            polygons = sockToWorldPolys(geometry, rad, plane, tol)
            polygons = [polygon for polygon in polygons if len(polygon) >= 3]
            if polygons:
                footprintGroups.append((index, polygons))
            else:
                rejectedIndices.append(index)
        worldPolys = [polygon for _, polygons in footprintGroups
                      for polygon in polygons]
        allPts = [pt for poly in worldPolys for pt in poly]
        if rad <= 0 and allPts:
            xs = [p[0] for p in allPts]
            ys = [p[1] for p in allPts]
            diag = math.hypot(max(xs)-min(xs), max(ys)-min(ys))
            rad = 0.03*diag if diag > 0 else 1.0
            footprintGroups = []
            rejectedIndices = list(unresolvedIndices)
            for index, geometry in resolvedPairs:
                polygons = sockToWorldPolys(geometry, rad, plane, tol)
                polygons = [polygon for polygon in polygons if len(polygon) >= 3]
                if polygons:
                    footprintGroups.append((index, polygons))
                else:
                    rejectedIndices.append(index)
            worldPolys = [polygon for _, polygons in footprintGroups
                          for polygon in polygons]
            allPts = [pt for poly in worldPolys for pt in poly]
        if not allPts:
            return (FlowMesh, FoodMesh,
                    ["No valid sock footprints (check the Socks type hint = Geometry)"],
                    Delivery, Deflection, Filtered)
        geometrySeconds = time.perf_counter() - tGeometry

        # --- grid, mask, solve ---
        tGrid = time.perf_counter()
        grid = buildGrid(allPts, flow, res, rad, domainPts)
        nx, ny, h = grid["nx"], grid["ny"], grid["h"]
        insideDomain = domainMask(grid, domainPts)
        solid = np.zeros((ny, nx), dtype=bool)
        representedObjects = 0
        representedInDomain = 0
        for _, polygons in footprintGroups:
            gridPolygons = worldPolysToGrid(polygons, grid)
            objectMask = rasterizePolygons(
                nx, ny, gridPolygons, ensure_each=True)
            if objectMask.any():
                representedObjects += 1
            if (objectMask & insideDomain).any():
                representedInDomain += 1
            solid |= objectMask
        # The supplied domain is the requested analysis site. Hidden objects
        # outside it must not alter the visible field.
        solid &= insideDomain
        displaySolid = solid | ~insideDomain
        band = bandAround(solid, self.BAND_W)
        flowSteps = min(self.STEPS, max(700, 8*max(nx, ny)))
        sestonSteps = max(600, min(2200, 8*max(nx, ny)))
        gridSeconds = time.perf_counter() - tGrid

        tFlow = time.perf_counter()
        ux, uy, speed = solveFlow(
            nx, ny, solid, self.U0, self.TAU, flowSteps,
            tol=1e-5, check=50)
        flowSeconds = time.perf_counter() - tFlow

        tSeston = time.perf_counter()
        seston = solveSeston(
            nx, ny, ux, uy, solid, band, upt, self.SESTON_D,
            steps=sestonSteps, tol=1e-6)
        sestonSeconds = time.perf_counter() - tSeston
        deficit = np.clip(1 - speed/self.U0, 0, 1)

        # field metrics for the legends (mean/max over fluid cells only)
        fl = ~solid & insideDomain
        anyFl = bool(fl.any())
        defMean = float(deficit[fl].mean()) if anyFl else 0.0
        defMax = float(deficit[fl].max()) if anyFl else 0.0
        removed = np.clip(1.0 - seston, 0, 1)
        remMean = float(removed[fl].mean()) if anyFl else 0.0
        remMax = float(removed[fl].max()) if anyFl else 0.0

        Wm, Hm = grid["Wm"], grid["Hm"]
        xh, yh = grid["xhat"], grid["yhat"]
        offV = (yh[0]*1.15*Hm, yh[1]*1.15*Hm)   # food panel offset, +yhat side

        tMeshes = time.perf_counter()
        FlowMesh = coloredMesh(deficit, cmapDeflect, displaySolid, grid, Z)
        FoodMesh = coloredMesh(seston, cmapFood, displaySolid, grid, Z)
        FoodMesh.Translate(rc.Geometry.Vector3d(offV[0], offV[1], 0.0))
        meshSeconds = time.perf_counter() - tMeshes

        # --- regions ---
        cellArea = h*h
        for ri, crv in enumerate(regionCrvs):
            rpoly = worldPolysToGrid([curvePoints(crv)], grid)
            rmask = rasterizePolygons(nx, ny, rpoly) & ~solid & insideDomain
            n = int(rmask.sum())
            if n == 0:
                Report.append("region %d: empty (no fluid cells inside)" % ri)
                Delivery.append(0.0)
                Deflection.append(0.0)
                Filtered.append(0.0)
                continue
            dlv = float((speed[rmask]/self.U0).mean())
            dfl = float(np.degrees(np.arctan2(np.abs(uy[rmask]),
                        np.maximum(ux[rmask], 1e-9))).mean())
            rem = float((1 - seston[rmask]).clip(0).mean())
            Delivery.append(round(dlv, 3))
            Deflection.append(round(dfl, 2))
            Filtered.append(round(rem, 4))
            Report.append("region %d | delivery %.2f | deflection %.1f deg | "
                          "food removed %.0f%% | patch area %.3f"
                          % (ri, dlv, dfl, rem*100, n*cellArea))
        if not regionCrvs:
            Report.append("No regions: add closed curves to score patches.")
        Report.append(
            "GEOMETRY | %d supplied | %d resolved | %d recognized | "
            "%d represented on grid | %d inside domain | "
            "%d footprints | %d solid cells | domain %s"
            % (len(suppliedItems), len(sockGeos), len(footprintGroups),
               representedObjects, representedInDomain,
               len(worldPolys), int(solid.sum()),
               "supplied" if domainPts else "automatic"))
        if rejectedIndices:
            Report.append(
                "GEOMETRY WARNING | rejected object indices: %s"
                % ", ".join(str(index) for index in rejectedIndices))
        if representedInDomain < representedObjects:
            Report.append(
                "DOMAIN WARNING | %d represented objects lie outside the "
                "supplied domain and are hidden from the panel."
                % (representedObjects-representedInDomain))

        # --- automatic viewport overlay: two side-by-side panels ---
        c00 = grid["cellToWorld"](-0.5, -0.5)
        c10 = grid["cellToWorld"](nx-0.5, -0.5)
        c11 = grid["cellToWorld"](nx-0.5, ny-0.5)
        c01 = grid["cellToWorld"](-0.5, ny-0.5)
        flowCorners = (c00, c10, c11, c01)
        foodCorners = tuple((c[0]+offV[0], c[1]+offV[1]) for c in flowCorners)
        csize = 0.03*math.hypot(Wm, Hm)

        bar1, frame1, texts1 = panelOverlay(
            flowCorners, cmapDeflect,
            "Flow slowdown  |  mean %.0f%%  max %.0f%%" % (defMean*100, defMax*100),
            [("1.0  stopped", 1.0), ("0.5  half speed", 0.5),
             ("0.0  full flow", 0.0)],
            Z, Wm, Hm, xh, yh, csize)
        bar2, frame2, texts2 = panelOverlay(
            foodCorners, cmapFood,
            "Food removed  |  mean %.0f%%  max %.0f%%" % (remMean*100, remMax*100),
            [("1.0  full food", 1.0), ("0.5", 0.5), ("0.0  starved", 0.0)],
            Z, Wm, Hm, xh, yh, csize)
        exactFlowFrame = domainFrame(domainPts, Z)
        exactFoodFrame = domainFrame(domainPts, Z, offV)
        if exactFlowFrame is not None:
            frame1 = exactFlowFrame
            frame2 = exactFoodFrame

        # incoming-flow arrow, upstream (-xh) of the flow panel, at mid-height
        midL = ((c00[0]+c01[0])*0.5, (c00[1]+c01[1])*0.5)
        aStart = rc.Geometry.Point3d(midL[0]-xh[0]*0.14*Wm, midL[1]-xh[1]*0.14*Wm, Z)
        aEnd = rc.Geometry.Point3d(midL[0]-xh[0]*0.02*Wm, midL[1]-xh[1]*0.02*Wm, Z)
        texts = texts1 + texts2 + [("FLOW", rc.Geometry.Point3d(
            aStart.X-xh[0]*0.04*Wm, aStart.Y-xh[1]*0.04*Wm, Z))]

        # clipping box over everything drawn (panels + legends + arrow)
        far = 0.15*Wm + 0.6*Hm
        pts = [rc.Geometry.Point3d(p[0], p[1], Z)
               for p in flowCorners + foodCorners]
        pts += [aStart, aEnd,
                rc.Geometry.Point3d(c10[0]+xh[0]*far, c10[1]+xh[1]*far, Z),
                rc.Geometry.Point3d(foodCorners[1][0]+xh[0]*far,
                                    foodCorners[1][1]+xh[1]*far, Z)]
        box = rc.Geometry.BoundingBox(pts)

        self._draw = dict(bars=[bar1, bar2], frames=[frame1, frame2],
                          arrow=rc.Geometry.Line(aStart, aEnd),
                          texts=texts, csize=csize, box=box)

        totalSeconds = time.perf_counter() - tAll
        Report.insert(0,
            "TIMING | grid %dx%d | geometry %.3fs | grid/mask %.3fs | "
            "flow %.3fs (%d max steps) | seston %.3fs (%d max steps) | "
            "Rhino meshes %.3fs | total %.3fs"
            % (nx, ny, geometrySeconds, gridSeconds, flowSeconds,
               flowSteps, sestonSeconds, sestonSteps,
               meshSeconds, totalSeconds))

        return (FlowMesh, FoodMesh, Report,
                Delivery, Deflection, Filtered)

    # ---- custom viewport preview (SDK-only) ----
    @property
    def ClippingBox(self):
        d = getattr(self, "_draw", None)
        if d:
            return d["box"]
        return rc.Geometry.BoundingBox.Empty

    def DrawViewportMeshes(self, args):
        """ Draw the two colour-scale legend bars (shaded via their vertex colors) """
        d = getattr(self, "_draw", None)
        if not d:
            return
        try:
            dp = args.Display
            for bar in d["bars"]:
                dp.DrawMeshFalseColors(bar)
        except Exception:
            pass

    def DrawViewportWires(self, args):
        """ Draw the panel frames, the incoming-flow arrow, and all text labels """
        d = getattr(self, "_draw", None)
        if not d:
            return
        try:
            dp = args.Display
            ink = System.Drawing.Color.FromArgb(35, 35, 35)
            arrows = d.get("arrows")
            if arrows is None and d.get("arrow") is not None:
                arrows = [d["arrow"]]
            for arrow in arrows or []:
                dp.DrawArrow(
                    arrow, System.Drawing.Color.FromArgb(0, 90, 200))
            for fr in d["frames"]:
                dp.DrawPolyline(fr, ink, 2)
            hgt = 0.45*d["csize"]
            for txt, pt in d["texts"]:
                plane = rc.Geometry.Plane(pt, rc.Geometry.Vector3d.ZAxis)
                dp.Draw3dText(txt, ink, plane, hgt)
        except Exception:
            pass
