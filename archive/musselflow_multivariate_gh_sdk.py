# r: numpy
"""
MusselFlow Multivariate Explorer
================================

Self-contained Rhino 8 Grasshopper Python 3 SDK-mode component. It accepts one
MusselFlow SiteData JSON document for a temporal analysis, or joins several
MusselFlow Site Field JSON documents by their stable grid-cell IDs for a spatial
analysis. It computes a Pearson correlation matrix and Principal Component
Analysis (PCA) with NumPy.

For the simplest temporal exploration, connect SiteData directly. For spatial
exploration, create several Site Field components with the same domain,
resolution, time, and depth but different variables, then connect all FieldData
outputs as one list here. Results describe associations between sampled fields;
they do not prove biological causation.

Inputs
------
run : item / bool
    True calculates the aligned multivariate analysis.
DataJson : list / str
    Either one SiteDataJson document for time-frame analysis, or two or more
    FieldDataJson documents from MusselFlow Site Field. Spatial fields should
    use the same domain, resolution, timeIndex, and depth.
components : item / int
    Requested PCA component count, clamped to the available rank. Three is a
    useful starting value.
standardize : item / bool
    True gives every variable unit variance before PCA. Recommended because
    oxygen, nutrients, temperature, and current use different physical units.
plotOrigin : item / Rhino.Geometry.Point3d
    Optional Rhino placement point for the correlation chart and PCA score plot.
cellSize : item / float
    Rhino size of one correlation-matrix cell. Defaults to 1.0.

Outputs — create ten ports once in this exact order
---------------------------------------------------
CorrelationMesh : item / Rhino.Geometry.Mesh
    Coloured Pearson correlation matrix: blue negative, white zero, red positive.
ScorePoints : list / Rhino.Geometry.Point3d
    PCA observations plotted as PC1 versus PC2 beside the correlation matrix.
Variables : list / str
    Variables retained after constant-field screening, in matrix column order.
CorrelationRows : list / str
    Human-readable correlation rows, each beginning with its variable name.
CorrelationMatrix : list / float
    Flattened row-major matrix. For p variables, every p values form one row.
Eigenvalues : list / float
    PCA eigenvalues in descending order.
ExplainedVariance : list / float
    Fraction of analyzed variance explained by each returned component.
Loadings : list / str
    One readable line per component listing each variable loading.
AnalysisDataJson : item / str
    Canonical JSON containing IDs, correlations, scores, loadings, and settings.
Report : list / str
    Alignment, exclusions, statistical limits, and interpretation guidance.

SDK setup
---------
Create a Rhino 8 Python 3 component, convert it with
``Convert To GH_ScriptInstance``, and replace its generated text with this file.
The RunScript annotations create and type-hint the six inputs. Add ten outputs
once in the exact order above. BeforeRunScript supplies names and hover help.
"""

import json
import math
import time

import Grasshopper
import Rhino
import System.Drawing
import numpy as np


COMPONENT_BUILD = "2026-08-08b"

COMPONENT_METADATA = {
    "name": "MusselFlow Multivariate Explorer",
    "nickname": "FieldPCA",
    "description": (
        "Analyze one temporal SiteData document or aligned spatial Site Field "
        "documents with Pearson correlation and PCA."),
}

INPUT_METADATA = (
    ("run", "run", "True calculates the correlation matrix and PCA."),
    ("DataJson", "Data", "One SiteDataJson, or aligned Site Field FieldDataJson documents."),
    ("components", "components", "Requested PCA component count; three recommended."),
    ("standardize", "standardize", "Standardize unlike physical units before PCA; recommended True."),
    ("plotOrigin", "plotOrigin", "Optional Rhino origin for the matrix and PCA score plot."),
    ("cellSize", "cellSize", "Rhino size of one correlation cell; defaults to 1.0."),
)

OUTPUT_METADATA = (
    ("CorrelationMesh", "CorrelationMesh", "Coloured Pearson correlation matrix."),
    ("ScorePoints", "ScorePoints", "PCA observations plotted as PC1 versus PC2."),
    ("Variables", "Variables", "Retained variables in matrix and loading order."),
    ("CorrelationRows", "CorrelationRows", "Readable correlation rows."),
    ("CorrelationMatrix", "CorrelationMatrix", "Flattened row-major Pearson matrix."),
    ("Eigenvalues", "Eigenvalues", "PCA eigenvalues in descending order."),
    ("ExplainedVariance", "ExplainedVariance", "Variance fraction explained by each PCA component."),
    ("Loadings", "Loadings", "Readable PCA loading lines."),
    ("AnalysisDataJson", "AnalysisData", "Canonical machine-readable correlation/PCA JSON."),
    ("Report", "Report", "Alignment, exclusions, cautions, and interpretation."),
)


def apply_component_metadata(component):
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "Correlation + PCA"
    for index, (name, nickname, description) in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description
        if index == 1:
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
    for index, (name, nickname, description) in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name = name
        parameter.NickName = nickname
        parameter.Description = description


def empty_outputs(status, message):
    return (None, [], [], [], [], [], [], [], "{}", [
        "MUSSELFLOW MULTIVARIATE | build %s | %s" % (COMPONENT_BUILD, status),
        message,
    ])


def canonical_json(value):
    return json.dumps(
        rounded(value), sort_keys=True, separators=(",", ":"),
        ensure_ascii=False, allow_nan=False)


def rounded(value, digits=6):
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, (float, np.floating)):
        return round(float(value), digits)
    if isinstance(value, np.ndarray):
        return rounded(value.tolist(), digits)
    if isinstance(value, (list, tuple)):
        return [rounded(item, digits) for item in value]
    if isinstance(value, dict):
        return {str(key): rounded(item, digits) for key, item in value.items()}
    return value


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def as_list(values):
    if values is None:
        return []
    if isinstance(values, (str, bytes)):
        return [values]
    try:
        return list(values)
    except TypeError:
        return [values]


def parse_field_document(source, index):
    try:
        document = json.loads(str(source))
    except Exception as exception:
        raise ValueError("fieldData[%d] is invalid JSON: %s" % (index, exception))
    if not isinstance(document, dict):
        raise ValueError("fieldData[%d] root is not an object." % index)
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("fieldData[%d] contains no records." % index)
    variable = str(document.get("variable") or "field_%d" % index)
    timestamp = str(document.get("time_utc") or "")
    units = str(document.get("units") or "unknown")
    values = {}
    points = {}
    for record in records:
        if not isinstance(record, dict):
            continue
        sample_id = record.get("sample_id")
        if sample_id is None:
            row = record.get("row")
            column = record.get("column")
            if row is None or column is None:
                continue
            sample_id = "%s:%s" % (row, column)
        sample_id = str(sample_id)
        value = finite_number(record.get("value"))
        if value is None:
            continue
        values[sample_id] = value
        point = record.get("rhino_point") or {}
        try:
            xyz = (
                float(point.get("x")),
                float(point.get("y")),
                float(point.get("z")))
            if all(math.isfinite(item) for item in xyz):
                points[sample_id] = xyz
        except Exception:
            pass
    if not values:
        raise ValueError("fieldData[%d] contains no finite values." % index)
    return {
        "variable": variable,
        "time": timestamp,
        "units": units,
        "values": values,
        "points": points,
        "grid": document.get("grid") or {},
        "depth_m": document.get("depth_m"),
    }


SITE_UNITS = {
    "chlorophyll_a_ug_l": "ug/L",
    "current_speed_m_s": "m/s",
    "dissolved_oxygen_mg_l": "mg/L",
    "eastward_current_m_s": "m/s",
    "nitrate_mmol_m3": "mmol/m3",
    "northward_current_m_s": "m/s",
    "phosphate_mmol_m3": "mmol/m3",
    "phytoplankton_carbon_mmol_m3": "mmol/m3",
    "salinity_psu": "PSU",
    "satellite_chlorophyll_a_ug_l": "ug/L",
    "temperature_c": "degC",
    "tsm_mg_l": "mg/L",
    "turbidity_fnu": "FNU",
}


def parse_json_document(source, index):
    try:
        document = json.loads(str(source))
    except Exception as exception:
        raise ValueError("data[%d] is invalid JSON: %s" % (index, exception))
    if not isinstance(document, dict):
        raise ValueError("data[%d] root is not an object." % index)
    return document


def temporal_site_matrix(document):
    frames = document.get("frames")
    if not isinstance(frames, list) or len(frames) < 3:
        raise ValueError("SiteData needs at least three valid time frames.")
    aggregates = []
    sample_ids = []
    for index, frame in enumerate(frames):
        aggregate = frame.get("aggregate") if isinstance(frame, dict) else None
        if not isinstance(aggregate, dict):
            raise ValueError("SiteData frame %d has no aggregate values." % index)
        aggregates.append(aggregate)
        sample_ids.append(str(frame.get("time_utc") or "frame_%d" % index))
    all_names = sorted(set().union(*(set(item) for item in aggregates)))
    labels = [
        name for name in all_names
        if all(finite_number(item.get(name)) is not None for item in aggregates)]
    missing = [name for name in all_names if name not in labels]
    if len(labels) < 2:
        raise ValueError(
            "SiteData has fewer than two variables complete across all frames.")
    matrix = np.asarray([
        [finite_number(aggregate[name]) for name in labels]
        for aggregate in aggregates], dtype=float)
    fields = []
    for column, label in enumerate(labels):
        fields.append({
            "variable": label,
            "time": "",
            "units": SITE_UNITS.get(label, "unknown"),
            "values": {
                sample_ids[row]: float(matrix[row, column])
                for row in range(len(sample_ids))},
            "points": {},
            "grid": {},
            "depth_m": None,
        })
    return fields, labels, sample_ids, matrix, missing


def unique_labels(fields):
    base_counts = {}
    for field in fields:
        base_counts[field["variable"]] = base_counts.get(field["variable"], 0)+1
    labels = []
    used = {}
    for field in fields:
        base = field["variable"]
        if base_counts[base] > 1:
            stamp = field["time"] or "undated"
            base = "%s@%s" % (base, stamp)
        used[base] = used.get(base, 0)+1
        label = base if used[base] == 1 else "%s#%d" % (base, used[base])
        labels.append(label)
    return labels


def aligned_matrix(fields):
    common = set(fields[0]["values"])
    union = set(common)
    for field in fields[1:]:
        keys = set(field["values"])
        common.intersection_update(keys)
        union.update(keys)
    sample_ids = sorted(common, key=sample_sort_key)
    if len(sample_ids) < 3:
        raise ValueError(
            "fewer than three common finite grid cells remain after alignment.")
    matrix = np.asarray([
        [field["values"][sample_id] for field in fields]
        for sample_id in sample_ids], dtype=float)
    return sample_ids, matrix, len(union)-len(common)


def sample_sort_key(value):
    try:
        row, column = str(value).split(":", 1)
        return int(row), int(column)
    except Exception:
        return str(value), ""


def remove_constant_fields(fields, labels, matrix):
    standard_deviation = np.std(matrix, axis=0, ddof=1)
    scale = np.maximum(np.max(np.abs(matrix), axis=0), 1.0)
    keep = standard_deviation > np.finfo(float).eps*scale*32.0
    excluded = [labels[index] for index in range(len(labels)) if not keep[index]]
    retained_fields = [fields[index] for index in range(len(fields)) if keep[index]]
    retained_labels = [labels[index] for index in range(len(labels)) if keep[index]]
    return retained_fields, retained_labels, matrix[:, keep], excluded


def calculate_pca(matrix, component_count, standardize):
    means = np.mean(matrix, axis=0)
    std = np.std(matrix, axis=0, ddof=1)
    centered = matrix-means
    transformed = centered/std if standardize else centered
    u, singular, vt = np.linalg.svd(transformed, full_matrices=False)
    rank_limit = min(matrix.shape[1], matrix.shape[0]-1)
    count = max(1, min(int(component_count), rank_limit, len(singular)))
    eigenvalues_all = (singular**2)/float(matrix.shape[0]-1)
    total = float(np.sum(eigenvalues_all))
    explained_all = (
        eigenvalues_all/total if total > np.finfo(float).eps
        else np.zeros_like(eigenvalues_all))
    return {
        "means": means,
        "standard_deviations": std,
        "scores": u[:, :count]*singular[:count],
        "loadings": vt[:count, :].T,
        "eigenvalues": eigenvalues_all[:count],
        "explained": explained_all[:count],
        "count": count,
    }


def diverging_color(value):
    t = max(-1.0, min(1.0, float(value)))
    white = (245, 245, 245)
    blue = (35, 90, 190)
    red = (210, 45, 40)
    endpoint = red if t >= 0.0 else blue
    amount = abs(t)
    rgb = [int(round(white[i]+(endpoint[i]-white[i])*amount)) for i in range(3)]
    return System.Drawing.Color.FromArgb(rgb[0], rgb[1], rgb[2])


def add_matrix_cell(mesh, origin, x, y, size, color):
    start = mesh.Vertices.Count
    z = origin.Z
    points = (
        Rhino.Geometry.Point3d(origin.X+x*size, origin.Y+y*size, z),
        Rhino.Geometry.Point3d(origin.X+(x+1)*size, origin.Y+y*size, z),
        Rhino.Geometry.Point3d(origin.X+(x+1)*size, origin.Y+(y+1)*size, z),
        Rhino.Geometry.Point3d(origin.X+x*size, origin.Y+(y+1)*size, z),
    )
    for point in points:
        mesh.Vertices.Add(point)
        mesh.VertexColors.Add(color.R, color.G, color.B)
    mesh.Faces.AddFace(int(start), int(start+1), int(start+2), int(start+3))


def matrix_mesh(correlation, origin, size):
    mesh = Rhino.Geometry.Mesh()
    count = correlation.shape[0]
    for row in range(count):
        for column in range(count):
            add_matrix_cell(
                mesh, origin, column, count-1-row, size,
                diverging_color(correlation[row, column]))
    mesh.Normals.ComputeNormals()
    mesh.Compact()
    return mesh


def score_points(scores, origin, matrix_width, cell_size):
    x = scores[:, 0]
    y = scores[:, 1] if scores.shape[1] > 1 else np.zeros(scores.shape[0])
    maximum = max(float(np.max(np.abs(x))), float(np.max(np.abs(y))), 1e-12)
    scale = max(2.0*cell_size, matrix_width*cell_size*0.45)/maximum
    offset_x = origin.X+(matrix_width+2.0)*cell_size
    offset_y = origin.Y+0.5*matrix_width*cell_size
    return [
        Rhino.Geometry.Point3d(
            offset_x+float(x[index])*scale,
            offset_y+float(y[index])*scale,
            origin.Z)
        for index in range(len(x))]


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            run: bool,
            data: list[str],
            components: int,
            standardize: bool,
            plotOrigin: Rhino.Geometry.Point3d,
            cellSize: float):
        """Join aligned physical fields and calculate correlation plus PCA."""
        started = time.perf_counter()
        if not run:
            return empty_outputs("WAITING", "Set run to True to calculate.")
        try:
            sources = [item for item in as_list(data)
                       if item is not None and str(item).strip()]
            if not sources:
                raise ValueError("connect SiteData or a list of Site Field FieldData JSON.")
            documents = [parse_json_document(source, index)
                         for index, source in enumerate(sources)]
            schemas = [str(document.get("schema") or "") for document in documents]
            if len(documents) == 1 and schemas[0].startswith("musselflow.site_data."):
                fields, labels, sample_ids, matrix, missing_fields = (
                    temporal_site_matrix(documents[0]))
                dropped_rows = 0
                analysis_scope = "temporal_site_frames"
            else:
                if any(not schema.startswith("musselflow.site_field.")
                       for schema in schemas):
                    raise ValueError(
                        "use one SiteData document, or only Site Field FieldData documents.")
                if len(sources) < 2:
                    raise ValueError(
                        "spatial mode needs at least two Site Field FieldData documents.")
                fields = [parse_field_document(source, index)
                          for index, source in enumerate(sources)]
                labels = unique_labels(fields)
                sample_ids, matrix, dropped_rows = aligned_matrix(fields)
                missing_fields = []
                analysis_scope = "spatial_grid_cell_association"
            fields, labels, matrix, excluded = remove_constant_fields(
                fields, labels, matrix)
            if len(labels) < 2:
                raise ValueError(
                    "fewer than two varying fields remain. The selected domain "
                    "may resolve to one repeated Copernicus source cell.")
            requested = 3 if components is None else int(components)
            requested = max(1, requested)
            use_standardization = True if standardize is None else bool(standardize)
            origin = Rhino.Geometry.Point3d.Origin
            try:
                candidate = Rhino.Geometry.Point3d(plotOrigin)
                if candidate.IsValid:
                    origin = candidate
            except Exception:
                pass
            size = finite_number(cellSize)
            size = 1.0 if size is None or size <= 0.0 else size
            correlation = np.corrcoef(matrix, rowvar=False)
            correlation = np.clip(correlation, -1.0, 1.0)
            pca = calculate_pca(matrix, requested, use_standardization)
        except Exception as exception:
            return empty_outputs("INVALID_INPUT", str(exception))

        correlation_mesh = matrix_mesh(correlation, origin, size)
        plotted_scores = score_points(
            pca["scores"], origin, len(labels), size)
        correlation_rows = [
            "%s | %s" % (
                labels[row], " | ".join(
                    "%+.4f" % correlation[row, column]
                    for column in range(len(labels))))
            for row in range(len(labels))]
        loading_lines = []
        for component in range(pca["count"]):
            loading_lines.append(
                "PC%d | %s" % (
                    component+1,
                    " | ".join(
                        "%s %+.4f" % (
                            labels[index], pca["loadings"][index, component])
                        for index in range(len(labels)))))

        original_points = {}
        for field in fields:
            original_points.update(field["points"])
        analysis_document = {
            "schema": "musselflow.multivariate.1.0.0",
            "build": COMPONENT_BUILD,
            "analysis_scope": analysis_scope,
            "variables": [
                {"name": labels[index], "units": fields[index]["units"],
                 "time_utc": fields[index]["time"],
                 "depth_m": fields[index]["depth_m"]}
                for index in range(len(labels))],
            "sample_ids": sample_ids,
            "sample_rhino_points": [
                original_points.get(sample_id) for sample_id in sample_ids],
            "standardized_pca": use_standardization,
            "correlation": correlation,
            "pca": {
                "eigenvalues": pca["eigenvalues"],
                "explained_variance": pca["explained"],
                "loadings": pca["loadings"],
                "scores": pca["scores"],
                "means": pca["means"],
                "standard_deviations": pca["standard_deviations"],
            },
            "excluded_constant_fields": excluded,
            "excluded_missing_fields": missing_fields,
            "dropped_unmatched_cells": dropped_rows,
        }
        report = [
            "MUSSELFLOW MULTIVARIATE | build %s | %d observations x %d variables | %.3f ms"
            % (COMPONENT_BUILD, len(sample_ids), len(labels),
               (time.perf_counter()-started)*1000.0),
            ("ROWS | %d time frames" % len(sample_ids)
             if analysis_scope == "temporal_site_frames"
             else "ALIGNMENT | %d common finite grid cells | %d unmatched cells excluded"
             % (len(sample_ids), dropped_rows)),
            "CORRELATION | Pearson coefficients from -1 to +1; diagonal is 1",
            "PCA | %d components | %s"
            % (pca["count"],
               "standardized to equal variance" if use_standardization
               else "covariance PCA in physical units"),
            "EXPLAINED | %s" % " | ".join(
                "PC%d %.2f%%" % (index+1, value*100.0)
                for index, value in enumerate(pca["explained"])),
        ]
        if excluded:
            report.append(
                "CONSTANT FIELD WARNING | excluded: %s" % ", ".join(excluded))
        if missing_fields:
            report.append(
                "MISSING FIELD WARNING | excluded because at least one frame is null: %s"
                % ", ".join(missing_fields))
        times = sorted(set(field["time"] for field in fields))
        depths = sorted(set(str(field["depth_m"]) for field in fields))
        if len(times) > 1:
            report.append(
                "TIME WARNING | inputs use different timestamps; spatial rows are "
                "aligned, but interpretation mixes environmental states.")
        if len(depths) > 1:
            report.append(
                "DEPTH WARNING | inputs use different depths; interpret associations carefully.")
        if analysis_scope == "temporal_site_frames":
            report.extend([
                "SAMPLE WARNING | fewer than 20 frames is exploratory only; fetch a longer period.",
                "TIME LIMIT | neighbouring frames are autocorrelated and not independent samples.",
            ])
        else:
            report.extend([
                "SCOPE | this run compares spatial grid cells, not an independent time series.",
                "RESOLUTION LIMIT | repeated regional source cells reduce effective sample size.",
            ])
        report.append(
            "CAUSALITY LIMIT | correlation and PCA reveal association, not cause and effect.")
        return (
            correlation_mesh,
            plotted_scores,
            labels,
            correlation_rows,
            rounded(correlation.reshape(-1).tolist()),
            rounded(pca["eigenvalues"].tolist()),
            rounded(pca["explained"].tolist()),
            loading_lines,
            canonical_json(analysis_document),
            report)

    def AfterRunScript(self):
        pass
