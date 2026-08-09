# MusselFlow — Current Grasshopper setup

There is one current editable ecological grammar and one current Grasshopper
component.

## Use these files

- `musselflow_component_gh_sdk.py` — paste this into one Rhino 8 Python 3
  Grasshopper component.
- `musselflow_grammar.json` — paste this complete JSON into a Grasshopper Panel
  connected to `caseJson`.
- `musselflow_case_schema.json` — machine-readable validation contract.
- `MUSSELFLOW_GRAMMAR_BOOK.md` — human-readable units, equations, validity
  envelopes, screening priors, and citations.

Keep these runtime sidecars beside the saved Grasshopper definition:

- `musselflow_case_core.py`
- `musselflow_ecogrammar_core.py`
- `musselflow_bio_optimizer_core.py`
- `musselflow_bio_optimizer_gh_sdk.py`

The last two sidecar names are retained temporarily for compatibility. They are
not alternative ecological grammars and should not be pasted as the current
component.

## Grasshopper component

1. Add a Rhino 8 Python 3 component.
2. Open its editor and choose `Convert To GH_ScriptInstance`.
3. Paste `musselflow_component_gh_sdk.py`.
4. Compile it so the typed signature synchronizes these inputs:

   - `run`
   - `obstacles`
   - `domain`
   - `probes`
   - `flowVectors`
   - `caseJson`

   `obstacles`, `probes`, and `flowVectors` use **List access**. They are not
   Item inputs and do not require grafting. Connect ordinary Grasshopper lists.

5. Give the component seven outputs in this order:

   - `Fitness`
   - `Feasible`
   - `Objectives`
   - `Constraints`
   - `Result`
   - `Status`
   - `Report`

6. Save the configured component as a Grasshopper User Object so its ports never
   need to be recreated.

## Scientific status

The component is a deterministic reduced-order comparative screening model.
It is not CFD, a calibrated carrying-capacity model, a legal assessment,
predicted growth, or verified nutrient/carbon removal.

All active ecological parameters are visible in `musselflow_grammar.json`.
The regression suite proves that it resolves to the complete 88-key numerical
configuration. The `evidence` records provide provenance but do not make a
screening prior site-calibrated.

## Archived material

Superseded split-version documents and the old Panel grammar are in `archive/`.
They exist only to explain or recover older Grasshopper definitions.
