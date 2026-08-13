# r: numpy
"""
MusselFlow PLS Field Lattice — 2.5D Environmental Relations
============================================================

Joins several MusselFlow Site Field outputs at one time and depth, fits a
supervised Partial Least Squares regression to one explicitly selected target,
then turns only statistically usable relationships into a vertical 2.5D Rhino
lattice. Every XY sample becomes a column of variable layers; local 3 x 3
correlations become angular links between the predictor and target layers.

The calculation runs before the drawing. Preview-grid duplicates are never
treated as independent evidence: PLS observations are collapsed to unique
Copernicus source cells of the target field, constant predictors are removed,
and the component stops when the effective sample is too small. Sentinel-3
OLCI surface fields and regional depth-aware fields may be compared, but their
different support, time and depth remain explicit in the result.

Name: MusselFlow PLS Field Lattice
Updated: 260813
Author: Felix Berger
Copyright: Apache License 2.0

Inputs:
run : item / bool
    True calculates PLS first and creates the 2.5D lattice afterwards.
FieldDataJson : list / str
    FieldDataJson outputs from MusselFlow Copernicus Field. Supply one document
    per variable, using the same Rhino domain, grid, timestamp and intended
    environmental state. Recorder placeholders, empty strings and malformed or
    unrelated JSON documents are skipped and identified in Report.
targetVariable : item / str
    Variable to predict, for example satellite_chlorophyll. PLS is supervised,
    so this target must be explicit and must vary across source cells.
maxComponents : item / int
    Maximum latent PLS components considered. Spatial block cross-validation
    chooses a value from 1 to this limit. Defaults to 3.
layerSpacing : item / float
    Rhino Z distance between variable layers. Defaults to one tenth of the XY
    field span when zero or missing.
thresholds : list / float
    Absolute local Pearson-r boundaries in descending order: high, medium and
    exploratory. Defaults to 0.85, 0.60 and 0.30. Values below the last boundary
    are omitted; tiers are effect-size labels, not significance tests.

Outputs — create ten ports once in this exact order:
PLSModelJson : item / str
    Model coefficients, VIP, spatial CV diagnostics, correlation matrix,
    source-cell support and local-link records.
Predictors : list / str
    Non-constant predictor variables retained by PLS.
Coefficients : list / float
    PLS coefficients in original physical units, parallel to Predictors.
VIP : list / float
    Variable Importance in Projection values, parallel to Predictors.
CorrelationMatrix : list / float
    Flattened Pearson matrix for target plus retained predictors.
CorrelationLinks : tree / Rhino.Geometry.Curve
    Three branches: {0} high, {1} medium, {2} exploratory local associations.
LatticePoints : tree / Rhino.Geometry.Point3d
    One branch per variable layer, created only after a valid PLS fit.
LatticeColors : tree / System.Drawing.Color
    Colours parallel to LatticePoints, normalized within each displayed layer.
LatticeLines : list / Rhino.Geometry.Curve
    Vertical lines connecting variable layers at each aligned XY sample.
Report : list / str
    Effective sample size, exclusions, CV diagnostics and scientific limits.

SDK setup:
Create a Rhino 8 Python 3 component, convert it with
``Convert To GH_ScriptInstance``, and replace its generated text with this file.
RunScript annotations create and type-hint the six inputs. Add ten outputs once
in the exact order above; BeforeRunScript applies names and hover descriptions.
"""

import json
import math
import time

import Grasshopper
import Rhino
import System.Drawing
import numpy as np


COMPONENT_BUILD = "2026-08-13b"
MIN_EFFECTIVE_OBSERVATIONS = 8
MIN_LOCAL_OBSERVATIONS = 6

COMPONENT_METADATA = {
    "name": "MusselFlow PLS Field Lattice",
    "nickname": "FieldPLS",
    "description": (
        "Fit guarded spatial PLS first, then draw its retained environmental "
        "relationships as a 2.5D variable lattice."),
}

INPUT_METADATA = (
    ("run", "run", "True fits PLS and then constructs the 2.5D lattice."),
    ("FieldDataJson", "Fields", "Aligned Site Field JSON documents. Empty, malformed and unrelated Recorder entries are skipped and reported."),
    ("targetVariable", "target", "Explicit response variable predicted by PLS, for example satellite_chlorophyll."),
    ("maxComponents", "maxComp", "Maximum PLS components; spatial block CV chooses the final count."),
    ("layerSpacing", "Z", "Z distance between variable layers; zero selects an automatic spacing."),
    ("thresholds", "tiers", "High, medium, exploratory absolute-r boundaries; default 0.85, 0.60, 0.30."),
)

OUTPUT_METADATA = (
    ("PLSModelJson", "PLSModel", "Canonical PLS model, validation and local-association JSON."),
    ("Predictors", "Predictors", "Predictors retained after constant/source-support screening."),
    ("Coefficients", "Coefficients", "PLS coefficients in original physical units."),
    ("VIP", "VIP", "Variable Importance in Projection values."),
    ("CorrelationMatrix", "Correlation", "Flattened Pearson matrix: target followed by predictors."),
    ("CorrelationLinks", "Links", "Tree branches {0}=high, {1}=medium, {2}=exploratory."),
    ("LatticePoints", "Points2.5D", "One point-tree branch per variable layer."),
    ("LatticeColors", "Colors2.5D", "Colours parallel to the 2.5D point tree."),
    ("LatticeLines", "Verticals", "Vertical curves connecting aligned variable layers."),
    ("Report", "Report", "PLS support, CV diagnostics, exclusions and scientific limits."),
)


def apply_component_metadata(component):
    if component is None:
        return
    component.Name = COMPONENT_METADATA["name"]
    component.NickName = COMPONENT_METADATA["nickname"]
    component.Description = COMPONENT_METADATA["description"]
    component.Message = "PLS -> 2.5D"
    for index, metadata in enumerate(INPUT_METADATA):
        if index >= component.Params.Input.Count:
            break
        parameter = component.Params.Input[index]
        parameter.Name, parameter.NickName, parameter.Description = metadata
        if index in (1, 5):
            parameter.Access = Grasshopper.Kernel.GH_ParamAccess.list
    for index, metadata in enumerate(OUTPUT_METADATA):
        if index >= component.Params.Output.Count:
            break
        parameter = component.Params.Output[index]
        parameter.Name, parameter.NickName, parameter.Description = metadata


def empty_tree():
    return Grasshopper.DataTree[object]()


def empty_outputs(status, message):
    return (
        "{}", [], [], [], [], empty_tree(), empty_tree(), empty_tree(), [],
        [
            "MUSSELFLOW PLS FIELD | build %s | %s" % (COMPONENT_BUILD, status),
            message,
        ])


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


def canonical_json(value):
    return json.dumps(
        rounded(value), sort_keys=True, separators=(",", ":"),
        ensure_ascii=False, allow_nan=False)


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
    if hasattr(values, "AllData"):
        try:
            return list(values.AllData())
        except Exception:
            pass
    try:
        return list(values)
    except TypeError:
        return [values]


def json_object(source, index):
    text = str(source).strip()
    if not text:
        raise ValueError("FieldDataJson[%d] is empty." % index)
    start = text.find("{")
    end = text.rfind("}")
    if start < 0 or end < start:
        raise ValueError("FieldDataJson[%d] does not contain a JSON object." % index)
    try:
        value = json.loads(text[start:end+1])
    except Exception as exception:
        raise ValueError("FieldDataJson[%d] is invalid JSON: %s" % (index, exception))
    if not isinstance(value, dict):
        raise ValueError("FieldDataJson[%d] root is not an object." % index)
    if not str(value.get("schema") or "").startswith("musselflow.site_field."):
        raise ValueError(
            "FieldDataJson[%d] is not a MusselFlow Site Field document." % index)
    return value


def source_key(record):
    latitude = finite_number(record.get("source_latitude"))
    longitude = finite_number(record.get("source_longitude"))
    if latitude is None or longitude is None:
        latitude = finite_number(record.get("latitude"))
        longitude = finite_number(record.get("longitude"))
    if latitude is None or longitude is None:
        return None
    return "%.7f,%.7f" % (latitude, longitude)


def parse_field(document, index):
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("FieldDataJson[%d] contains no records." % index)
    variable = str(document.get("variable") or "field_%d" % index).strip()
    parsed = {}
    for record in records:
        if not isinstance(record, dict):
            continue
        value = finite_number(record.get("value"))
        row = finite_number(record.get("row"))
        column = finite_number(record.get("column"))
        sample_id = record.get("sample_id")
        if sample_id is None and row is not None and column is not None:
            sample_id = "%d:%d" % (int(row), int(column))
        if sample_id is None or value is None:
            continue
        point_data = record.get("rhino_point") or {}
        point = None
        try:
            xyz = [float(point_data.get(axis)) for axis in ("x", "y", "z")]
            if all(math.isfinite(number) for number in xyz):
                point = tuple(xyz)
        except Exception:
            point = None
        parsed[str(sample_id)] = {
            "value": value,
            "point": point,
            "row": None if row is None else int(row),
            "column": None if column is None else int(column),
            "source_key": source_key(record),
            "source_latitude": finite_number(record.get("source_latitude")),
            "source_longitude": finite_number(record.get("source_longitude")),
        }
    if not parsed:
        raise ValueError("FieldDataJson[%d] contains no finite records." % index)
    return {
        "variable": variable,
        "units": str(document.get("units") or "unknown"),
        "time": str(document.get("time_utc") or ""),
        "depth_m": finite_number(document.get("depth_m")),
        "nominal_resolution_m": finite_number(document.get("nominal_resolution_m")),
        "source": str(document.get("data_source") or ""),
        "grid": document.get("grid") or {},
        "records": parsed,
    }


def valid_field_inputs(values):
    """Return parsed Site Field documents and rejection records.

    Grasshopper Recorders often preserve empty placeholders such as an empty
    JSON object or unrelated status strings. Those are not statistical
    observations and must not invalidate otherwise usable fields. Original
    input indices are retained so Report identifies every ignored entry.
    """
    sources = as_list(values)
    accepted = []
    skipped = []
    for source_index, source in enumerate(sources):
        if source is None or not str(source).strip():
            skipped.append({
                "index": source_index,
                "reason": "empty input",
            })
            continue
        try:
            document = json_object(source, source_index)
            field = parse_field(document, source_index)
        except Exception as exception:
            skipped.append({
                "index": source_index,
                "reason": str(exception),
            })
            continue
        accepted.append({
            "index": source_index,
            "document": document,
            "field": field,
        })
    return sources, accepted, skipped

def unique_field_names(fields):
    counts = {}
    for field in fields:
        counts[field["variable"]] = counts.get(field["variable"], 0)+1
    labels = []
    used = {}
    for field in fields:
        base = field["variable"]
        if counts[base] > 1:
            base = "%s@%s" % (base, field["time"] or "undated")
        used[base] = used.get(base, 0)+1
        labels.append(base if used[base] == 1 else "%s#%d" % (base, used[base]))
    return labels


def normalized_name(value):
    return "".join(character for character in str(value).lower() if character.isalnum())


def target_index(labels, requested):
    target = normalized_name(requested)
    if not target:
        raise ValueError(
            "targetVariable is required. PLS cannot run without a declared response.")
    exact = [index for index, label in enumerate(labels)
             if normalized_name(label) == target]
    if len(exact) == 1:
        return exact[0]
    prefix = [index for index, label in enumerate(labels)
              if normalized_name(label).startswith(target)]
    if len(prefix) == 1:
        return prefix[0]
    raise ValueError(
        "targetVariable '%s' does not uniquely match: %s"
        % (requested, ", ".join(labels)))


def align_fields(fields):
    common = set(fields[0]["records"])
    union = set(common)
    for field in fields[1:]:
        keys = set(field["records"])
        common.intersection_update(keys)
        union.update(keys)
    sample_ids = sorted(common, key=sample_sort_key)
    if len(sample_ids) < MIN_EFFECTIVE_OBSERVATIONS:
        raise ValueError(
            "only %d common display samples remain; at least %d are required."
            % (len(sample_ids), MIN_EFFECTIVE_OBSERVATIONS))
    matrix = np.asarray([
        [field["records"][sample_id]["value"] for field in fields]
        for sample_id in sample_ids], dtype=float)
    return sample_ids, matrix, len(union)-len(common)


def sample_sort_key(value):
    try:
        row, column = str(value).split(":", 1)
        return int(row), int(column)
    except Exception:
        return str(value), ""


def aggregate_target_cells(fields, sample_ids, matrix, target):
    target_records = fields[target]["records"]
    groups = {}
    for row_index, sample_id in enumerate(sample_ids):
        record = target_records[sample_id]
        key = record["source_key"]
        if key is None:
            continue
        groups.setdefault(key, []).append(row_index)
    if len(groups) < MIN_EFFECTIVE_OBSERVATIONS:
        raise ValueError(
            "%d preview samples collapse to only %d independent target source "
            "cell(s); PLS needs at least %d. Enlarge the georeferenced sampling "
            "domain or add genuine time observations. Preview resolution cannot "
            "create statistical evidence."
            % (len(sample_ids), len(groups), MIN_EFFECTIVE_OBSERVATIONS))
    keys = sorted(groups)
    effective = []
    coordinates = []
    membership = []
    for key in keys:
        indices = groups[key]
        effective.append(np.mean(matrix[indices, :], axis=0))
        first = target_records[sample_ids[indices[0]]]
        latitude = first["source_latitude"]
        longitude = first["source_longitude"]
        coordinates.append((latitude, longitude))
        membership.append([sample_ids[index] for index in indices])
    return keys, np.asarray(effective, dtype=float), coordinates, membership


def varying_columns(matrix, target):
    deviations = np.std(matrix, axis=0, ddof=1)
    scale = np.maximum(np.max(np.abs(matrix), axis=0), 1.0)
    varying = deviations > np.finfo(float).eps*scale*64.0
    if not varying[target]:
        raise ValueError("targetVariable is constant across independent source cells.")
    predictors = [index for index in range(matrix.shape[1])
                  if index != target and varying[index]]
    excluded = [index for index in range(matrix.shape[1])
                if index != target and not varying[index]]
    if not predictors:
        raise ValueError("all predictors are constant across independent source cells.")
    return predictors, excluded


def fit_pls1(x, y, component_count, response_scale=1.0):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    x_mean = np.mean(x, axis=0)
    y_mean = float(np.mean(y))
    x_scale = np.std(x, axis=0, ddof=1)
    y_scale = float(np.std(y, ddof=1)) if response_scale else 1.0
    if np.any(x_scale <= np.finfo(float).eps) or y_scale <= np.finfo(float).eps:
        raise ValueError("PLS received a constant training field.")
    x_residual = (x-x_mean)/x_scale
    y_residual = (y-y_mean)/y_scale
    limit = min(int(component_count), x.shape[1], x.shape[0]-1)
    weights = []
    loadings = []
    y_loadings = []
    scores = []
    for unused in range(limit):
        weight = np.dot(x_residual.T, y_residual)
        norm = float(np.linalg.norm(weight))
        if norm <= np.finfo(float).eps:
            break
        weight = weight/norm
        score = np.dot(x_residual, weight)
        denominator = float(np.dot(score, score))
        if denominator <= np.finfo(float).eps:
            break
        loading = np.dot(x_residual.T, score)/denominator
        y_loading = float(np.dot(y_residual, score)/denominator)
        x_residual = x_residual-np.outer(score, loading)
        y_residual = y_residual-score*y_loading
        weights.append(weight)
        loadings.append(loading)
        y_loadings.append(y_loading)
        scores.append(score)
    if not weights:
        raise ValueError("PLS could not extract a non-zero latent component.")
    w = np.column_stack(weights)
    p = np.column_stack(loadings)
    q = np.asarray(y_loadings, dtype=float)
    t = np.column_stack(scores)
    rotations = np.dot(w, np.linalg.pinv(np.dot(p.T, w)))
    coefficient_standardized = np.dot(rotations, q)
    coefficients = coefficient_standardized*y_scale/x_scale
    intercept = y_mean-float(np.dot(x_mean, coefficients))
    prediction = np.dot(x, coefficients)+intercept
    explained_y = (q**2)*np.sum(t**2, axis=0)
    total_explained = float(np.sum(explained_y))
    if total_explained <= np.finfo(float).eps:
        vip = np.zeros(x.shape[1])
    else:
        normalized_weights = w**2/np.sum(w**2, axis=0)
        vip = np.sqrt(
            x.shape[1]*np.dot(normalized_weights, explained_y)/total_explained)
    return {
        "components": w.shape[1],
        "coefficients": coefficients,
        "intercept": intercept,
        "prediction": prediction,
        "weights": w,
        "loadings": p,
        "scores": t,
        "y_loadings": q,
        "vip": vip,
        "x_mean": x_mean,
        "x_scale": x_scale,
        "y_mean": y_mean,
        "y_scale": y_scale,
    }


def predict_pls(model, x):
    return np.dot(np.asarray(x, dtype=float), model["coefficients"])+model["intercept"]


def spatial_folds(coordinates, desired=4):
    coordinates = np.asarray(coordinates, dtype=float)
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("target source coordinates are incomplete.")
    latitude_range = float(np.ptp(coordinates[:, 0]))
    longitude_range = float(np.ptp(coordinates[:, 1]))
    axis = 1 if longitude_range >= latitude_range else 0
    order = np.argsort(coordinates[:, axis], kind="mergesort")
    count = min(int(desired), max(3, len(order)//3))
    count = min(count, len(order))
    folds = []
    for indices in np.array_split(order, count):
        if len(indices):
            folds.append(np.asarray(indices, dtype=int))
    if len(folds) < 3:
        raise ValueError("fewer than three spatial folds can be formed.")
    return folds


def cross_validate_components(x, y, coordinates, maximum):
    folds = spatial_folds(coordinates)
    maximum = min(int(maximum), x.shape[1], x.shape[0]-2)
    if maximum < 1:
        raise ValueError("no valid PLS component count remains.")
    rows = []
    for count in range(1, maximum+1):
        fold_rmse = []
        predictions = np.full(len(y), np.nan)
        valid = True
        for test in folds:
            train_mask = np.ones(len(y), dtype=bool)
            train_mask[test] = False
            if np.sum(train_mask) <= count+1:
                valid = False
                break
            try:
                model = fit_pls1(
                    x[train_mask], y[train_mask], count,
                    response_scale=0.0)
                predicted = predict_pls(model, x[test])
            except Exception:
                valid = False
                break
            predictions[test] = predicted
            fold_rmse.append(float(np.sqrt(np.mean((predicted-y[test])**2))))
        if valid and np.all(np.isfinite(predictions)):
            rows.append({
                "components": count,
                "mean_rmse": float(np.mean(fold_rmse)),
                "se_rmse": (float(np.std(fold_rmse, ddof=1)/math.sqrt(len(fold_rmse)))
                            if len(fold_rmse) > 1 else 0.0),
                "fold_rmse": fold_rmse,
                "predictions": predictions,
            })
    if not rows:
        raise ValueError(
            "spatial block cross-validation failed; increase independent spatial support.")
    best = min(rows, key=lambda row: row["mean_rmse"])
    limit = best["mean_rmse"]+best["se_rmse"]
    selected = min(
        (row for row in rows if row["mean_rmse"] <= limit),
        key=lambda row: row["components"])
    total = float(np.sum((y-np.mean(y))**2))
    press = float(np.sum((selected["predictions"]-y)**2))
    q2 = None if total <= np.finfo(float).eps else 1.0-press/total
    return selected, rows, q2, len(folds)


def parse_thresholds(values):
    numbers = [finite_number(value) for value in as_list(values)]
    numbers = [number for number in numbers if number is not None]
    if not numbers:
        numbers = [0.85, 0.60, 0.30]
    if len(numbers) != 3:
        raise ValueError("thresholds must contain exactly three values.")
    numbers = sorted([max(0.0, min(1.0, number)) for number in numbers], reverse=True)
    if not (numbers[0] > numbers[1] > numbers[2] > 0.0):
        raise ValueError("thresholds must be three distinct descending values above zero.")
    return numbers


def tier_index(absolute_r, thresholds):
    if absolute_r >= thresholds[0]:
        return 0
    if absolute_r >= thresholds[1]:
        return 1
    if absolute_r >= thresholds[2]:
        return 2
    return None


def point_for(field, sample_id, z):
    point = field["records"][sample_id]["point"]
    if point is None:
        return None
    return Rhino.Geometry.Point3d(point[0], point[1], z)


def estimate_spacing(fields, sample_ids):
    points = [fields[0]["records"][sample_id]["point"] for sample_id in sample_ids]
    points = [point for point in points if point is not None]
    if not points:
        return 1.0, 1.0, 1.0
    xs = sorted(set(round(point[0], 9) for point in points))
    ys = sorted(set(round(point[1], 9) for point in points))
    dx = min((xs[index+1]-xs[index] for index in range(len(xs)-1)), default=1.0)
    dy = min((ys[index+1]-ys[index] for index in range(len(ys)-1)), default=1.0)
    span = max(max(xs)-min(xs), max(ys)-min(ys), dx, dy, 1.0)
    return abs(dx), abs(dy), span


def layer_color(value, minimum, maximum, layer):
    amount = 0.5 if maximum <= minimum else (value-minimum)/(maximum-minimum)
    amount = max(0.0, min(1.0, amount))
    hues = (0.58, 0.48, 0.15, 0.02, 0.78, 0.90, 0.32, 0.66)
    hue = (hues[layer % len(hues)]+0.13*(amount-0.5)) % 1.0
    return Rhino.Display.ColorHSL(hue, 0.82, 0.35+0.30*amount).ToArgbColor()


def data_tree(branches):
    tree = Grasshopper.DataTree[object]()
    for branch_index, values in enumerate(branches):
        path = Grasshopper.Kernel.Data.GH_Path(branch_index)
        for value in values:
            tree.Add(value, path)
    return tree


def angular_link(start, end, dx, dy, sign):
    direction = 1.0 if sign >= 0.0 else -1.0
    z_one = start.Z+(end.Z-start.Z)*0.34
    z_two = start.Z+(end.Z-start.Z)*0.68
    points = [
        start,
        Rhino.Geometry.Point3d(start.X+direction*0.28*dx, start.Y, z_one),
        Rhino.Geometry.Point3d(start.X+direction*0.28*dx, start.Y+0.28*dy, z_two),
        end,
    ]
    return Rhino.Geometry.PolylineCurve(points)


def local_links(fields, labels, sample_ids, target, predictor_indices, z_values,
                thresholds, dx, dy, allowed_predictors):
    lookup = {}
    for sample_id in sample_ids:
        record = fields[target]["records"][sample_id]
        if record["row"] is not None and record["column"] is not None:
            lookup[(record["row"], record["column"])] = sample_id
    branches = [[], [], []]
    records = []
    if not lookup:
        return branches, records
    rows = [key[0] for key in lookup]
    columns = [key[1] for key in lookup]
    for row in range(min(rows)+1, max(rows)):
        for column in range(min(columns)+1, max(columns)):
            center_id = lookup.get((row, column))
            if center_id is None:
                continue
            neighborhood = [lookup.get((r, c))
                            for r in range(row-1, row+2)
                            for c in range(column-1, column+2)]
            neighborhood = [sample_id for sample_id in neighborhood if sample_id]
            for predictor in predictor_indices:
                if predictor not in allowed_predictors:
                    continue
                pairs = {}
                for sample_id in neighborhood:
                    target_record = fields[target]["records"][sample_id]
                    predictor_record = fields[predictor]["records"][sample_id]
                    pair_key = (target_record["source_key"], predictor_record["source_key"])
                    if None in pair_key:
                        continue
                    pairs.setdefault(pair_key, []).append((
                        target_record["value"], predictor_record["value"]))
                observations = [np.mean(values, axis=0) for values in pairs.values()]
                if len(observations) < MIN_LOCAL_OBSERVATIONS:
                    continue
                array = np.asarray(observations, dtype=float)
                if np.std(array[:, 0], ddof=1) <= np.finfo(float).eps:
                    continue
                if np.std(array[:, 1], ddof=1) <= np.finfo(float).eps:
                    continue
                correlation = float(np.corrcoef(array[:, 0], array[:, 1])[0, 1])
                tier = tier_index(abs(correlation), thresholds)
                if tier is None:
                    continue
                start = point_for(fields[predictor], center_id, z_values[predictor])
                end = point_for(fields[target], center_id, z_values[target])
                if start is None or end is None:
                    continue
                curve = angular_link(start, end, dx, dy, correlation)
                branches[tier].append(curve)
                records.append({
                    "sample_id": center_id,
                    "row": row,
                    "column": column,
                    "predictor": labels[predictor],
                    "target": labels[target],
                    "pearson_r": correlation,
                    "tier": ("high", "medium", "exploratory")[tier],
                    "independent_pairs": len(observations),
                })
    return branches, records


def build_lattice(fields, labels, sample_ids, target, predictor_indices,
                  spacing, thresholds, allowed_predictors):
    displayed_predictors = [index for index in predictor_indices
                            if index in allowed_predictors]
    displayed_indices = [target]+displayed_predictors
    dx, dy, span = estimate_spacing(fields, sample_ids)
    z_step = spacing if spacing is not None and spacing > 0.0 else 0.10*span
    z_base = 0.0
    for sample_id in sample_ids:
        point = fields[target]["records"][sample_id]["point"]
        if point is not None:
            z_base = point[2]
            break
    z_values = {field_index: z_base+layer*z_step
                for layer, field_index in enumerate(displayed_indices)}
    point_branches = []
    color_branches = []
    for layer, field_index in enumerate(displayed_indices):
        values = [fields[field_index]["records"][sample_id]["value"]
                  for sample_id in sample_ids]
        minimum = min(values)
        maximum = max(values)
        points = []
        colors = []
        for sample_id, value in zip(sample_ids, values):
            point = point_for(fields[field_index], sample_id, z_values[field_index])
            if point is None:
                continue
            points.append(point)
            colors.append(layer_color(value, minimum, maximum, layer))
        point_branches.append(points)
        color_branches.append(colors)
    verticals = []
    if len(displayed_indices) > 1:
        low = z_values[displayed_indices[0]]
        high = z_values[displayed_indices[-1]]
        for sample_id in sample_ids:
            point = point_for(fields[target], sample_id, low)
            if point is not None:
                verticals.append(Rhino.Geometry.LineCurve(
                    point, Rhino.Geometry.Point3d(point.X, point.Y, high)))
    link_branches, link_records = local_links(
        fields, labels, sample_ids, target, predictor_indices, z_values,
        thresholds, dx, dy, allowed_predictors)
    return (
        data_tree(point_branches), data_tree(color_branches), verticals,
        data_tree(link_branches), link_records, z_step, displayed_indices)


class Script_Instance(Grasshopper.Kernel.GH_ScriptInstance):

    def BeforeRunScript(self):
        apply_component_metadata(getattr(self, "Component", None))

    def RunScript(
            self,
            run: bool,
            FieldDataJson: list[str],
            targetVariable: str,
            maxComponents: int,
            layerSpacing: float,
            thresholds: list[float]):
        """Fit spatial PLS, then display retained relations as a 2.5D lattice."""
        started = time.perf_counter()
        if not run:
            return empty_outputs("WAITING", "Set run to True to calculate.")
        try:
            sources, accepted_inputs, skipped_inputs = valid_field_inputs(
                FieldDataJson)
            fields = [item["field"] for item in accepted_inputs]
            if len(fields) < 3:
                reasons = " | ".join(
                    "[%d] %s" % (item["index"], item["reason"])
                    for item in skipped_inputs)
                raise ValueError(
                    "need at least three valid MusselFlow Site Field documents "
                    "(one target and two predictors); received %d input(s), "
                    "accepted %d, skipped %d.%s"
                    % (len(sources), len(fields), len(skipped_inputs),
                       "" if not reasons else " SKIPPED | "+reasons))
            labels = unique_field_names(fields)
            target = target_index(labels, targetVariable)
            sample_ids, display_matrix, unmatched = align_fields(fields)
            source_keys, effective_matrix, source_coordinates, membership = (
                aggregate_target_cells(fields, sample_ids, display_matrix, target))
            predictor_indices, excluded_indices = varying_columns(
                effective_matrix, target)
            x = effective_matrix[:, predictor_indices]
            y = effective_matrix[:, target]
            requested_components = 3 if maxComponents is None else int(maxComponents)
            requested_components = max(1, requested_components)
            selected_cv, cv_rows, q2, fold_count = cross_validate_components(
                x, y, source_coordinates, requested_components)
            model = fit_pls1(
                x, y, selected_cv["components"], response_scale=0.0)
            retained_predictor_indices = {
                predictor_indices[index]
                for index, vip in enumerate(model["vip"])
                if vip >= 1.0
            }
            if not retained_predictor_indices:
                retained_predictor_indices = {
                    predictor_indices[int(np.argmax(model["vip"]))]
                }
            correlation_order = [target]+predictor_indices
            correlation = np.corrcoef(
                effective_matrix[:, correlation_order], rowvar=False)
            correlation = np.clip(correlation, -1.0, 1.0)
            tier_values = parse_thresholds(thresholds)
            spacing = finite_number(layerSpacing)
            lattice_points, lattice_colors, verticals, links, link_records, z_step, displayed = (
                build_lattice(
                    fields, labels, sample_ids, target, predictor_indices,
                    spacing, tier_values, retained_predictor_indices))
        except Exception as exception:
            return empty_outputs("INVALID_INPUT", str(exception))

        residual = y-model["prediction"]
        rmse_fit = float(np.sqrt(np.mean(residual**2)))
        total = float(np.sum((y-np.mean(y))**2))
        r2_fit = None if total <= np.finfo(float).eps else 1.0-float(np.sum(residual**2))/total
        predictor_labels = [labels[index] for index in predictor_indices]
        excluded_labels = [labels[index] for index in excluded_indices]
        target_field = fields[target]
        model_document = {
            "schema": "musselflow.pls_field_lattice.1.0.0",
            "build": COMPONENT_BUILD,
            "calculation_order": ["aligned_source_observations", "spatial_pls", "local_3x3_correlations", "2.5d_lattice"],
            "target": {
                "name": labels[target],
                "units": target_field["units"],
                "source": target_field["source"],
                "nominal_resolution_m": target_field["nominal_resolution_m"],
            },
            "predictors": [
                {
                    "name": labels[index],
                    "units": fields[index]["units"],
                    "coefficient": model["coefficients"][position],
                    "vip": model["vip"][position],
                    "source": fields[index]["source"],
                    "nominal_resolution_m": fields[index]["nominal_resolution_m"],
                }
                for position, index in enumerate(predictor_indices)],
            "pls": {
                "algorithm": "PLS1 NIPALS, X centered/unit-variance scaled; Y centered in physical target units",
                "components": model["components"],
                "intercept": model["intercept"],
                "fit_rmse": rmse_fit,
                "fit_r2": r2_fit,
                "spatial_cv_q2": q2,
                "spatial_cv_rmse": selected_cv["mean_rmse"],
                "spatial_cv_folds": fold_count,
                "component_search": [
                    {
                        "components": row["components"],
                        "mean_rmse": row["mean_rmse"],
                        "se_rmse": row["se_rmse"],
                        "fold_rmse": row["fold_rmse"],
                    }
                    for row in cv_rows],
                "selection_rule": "smallest component count within one SE of minimum spatial-CV RMSE",
            },
            "correlation": {
                "variables": [labels[index] for index in correlation_order],
                "pearson_matrix": correlation,
                "local_thresholds_absolute_r": {
                    "high": tier_values[0],
                    "medium": tier_values[1],
                    "exploratory": tier_values[2],
                },
                "local_3x3_links": link_records,
                "local_minimum_independent_pairs": MIN_LOCAL_OBSERVATIONS,
            },
            "support": {
                "display_samples": len(sample_ids),
                "effective_target_source_cells": len(source_keys),
                "target_source_cells": source_keys,
                "display_membership_by_target_cell": membership,
                "unmatched_display_samples": unmatched,
                "excluded_constant_predictors": excluded_labels,
            },
            "lattice": {
                "variable_order_bottom_to_top": [labels[index] for index in displayed],
                "selection": "target plus PLS predictors with VIP >= 1.0; strongest VIP retained if none reach 1.0",
                "layer_spacing_rhino_units": z_step,
                "display_only": True,
            },
            "input_filter": {
                "received": len(sources),
                "accepted": len(accepted_inputs),
                "accepted_original_indices": [
                    item["index"] for item in accepted_inputs],
                "skipped": skipped_inputs,
            },
            "scientific_status": "EXPLORATORY_SPATIAL_ASSOCIATION_NOT_CAUSATION_OR_SUPER_RESOLUTION",
        }
        link_counts = [len(link_records_for_branch)
                       for link_records_for_branch in links.Branches]
        while len(link_counts) < 3:
            link_counts.append(0)
        times = sorted(set(field["time"] for field in fields))
        depths = sorted(set(str(field["depth_m"]) for field in fields))
        report = [
            "MUSSELFLOW PLS FIELD | build %s | %d display samples -> %d independent target cells | %.3f ms"
            % (COMPONENT_BUILD, len(sample_ids), len(source_keys),
               (time.perf_counter()-started)*1000.0),
            "TARGET | %s | %s | source support %s m"
            % (labels[target], target_field["units"],
               "unknown" if target_field["nominal_resolution_m"] is None
               else "%.6g" % target_field["nominal_resolution_m"]),
            "PLS | %d component(s) selected by %d-fold spatial block CV | RMSECV %.6f | Q2 %s"
            % (model["components"], fold_count, selected_cv["mean_rmse"],
               "undefined" if q2 is None else "%.6f" % q2),
            "PREDICTORS | %s" % " | ".join(
                "%s coef %+.6f VIP %.6f"
                % (label, model["coefficients"][index], model["vip"][index])
                for index, label in enumerate(predictor_labels)),
            "LOCAL LINKS | high %d | medium %d | exploratory %d | thresholds %.2f / %.2f / %.2f"
            % (link_counts[0], link_counts[1], link_counts[2],
               tier_values[0], tier_values[1], tier_values[2]),
            "ORDER | PLS and correlation screening completed before 2.5D lattice construction.",
            "LATTICE FILTER | target plus PLS-retained predictors: %s"
            % " | ".join(labels[index] for index in displayed),
            "SOURCE-CELL GUARD | preview duplicates were collapsed by the target's Copernicus source cell.",
            "INPUT FILTER | %d received | %d valid Site Field document(s) | %d skipped"
            % (len(sources), len(accepted_inputs), len(skipped_inputs)),
        ]
        report.extend(
            "SKIPPED INPUT [%d] | %s" % (item["index"], item["reason"])
            for item in skipped_inputs)
        if excluded_labels:
            report.append("CONSTANT PREDICTORS EXCLUDED | "+" | ".join(excluded_labels))
        if len(times) > 1:
            report.append(
                "TIME WARNING | inputs use different timestamps: %s" % " | ".join(times))
        if len(depths) > 1:
            report.append(
                "DEPTH WARNING | inputs use different depths: %s" % " | ".join(depths))
        report.extend([
            "SENTINEL-3 LIMIT | OLCI ocean-colour fields describe the optically observed surface at native 300 m support and require cloud/quality screening.",
            "RESOLUTION LIMIT | the lattice may draw densely, but it does not create 10 m measurements or independent observations.",
            "STATISTICAL LIMIT | local tiers are absolute Pearson effect sizes, not p-values; many spatial comparisons remain exploratory.",
            "CAUSALITY LIMIT | PLS and correlation model association, not ecological causation.",
        ])
        return (
            canonical_json(model_document),
            predictor_labels,
            rounded(model["coefficients"]),
            rounded(model["vip"]),
            rounded(correlation.reshape(-1)),
            links,
            lattice_points,
            lattice_colors,
            verticals,
            report)

    def AfterRunScript(self):
        pass
