# Auto Bonded Contacts for Mecway

IronPython script that finds mating free faces between solid and shell components and writes bonded contacts straight into the Mecway `.liml` file.

Mecway's scripting API can read the mesh, but it cannot create contact pairs. This script therefore saves the model, keeps a `.bak` copy, and inserts `<faceselection>` and `<contact2 type="bondedce">` elements into the XML. After the run you must reload the file. Saving again before the reload overwrites the contacts.

## What you need

- Mecway with the IronPython scripting host (tested against `.liml` version 36).
- A saved model. Coordinates in the file are metres. The tolerance you type is millimetres and is converted internally.
- At least two components that contain solid or shell elements. Beams are ignored.
- Components must already exist. The script does not split a mesh into components.

Supported shapes:

| Role | Shapes | Faces used |
| --- | --- | --- |
| Solid | `tet4`, `tet10`, `pyr5`, `pyr13`, `wedge6`, `wedge15`, `penta6`, `penta15`, `hex8`, `hex20` | All element faces. Only faces whose node set appears once inside the component (free faces) are kept. |
| Shell | `tri3`, `tri6`, `quad4`, `quad8` | Surface faces 5 and 6 only. The side whose normal points at the partner is stored. Edges are not contacts. |
| Ignored | `line*`, `beam`, `truss`, `spring`, `damper`, anything else | Not used. |

## Install

1. Download `AutoBondedContacts.py` from the [latest release](https://github.com/larsmei/mecway-auto-bonded-contacts/releases/latest).
2. In Mecway choose **Tools → Scripts → New** and paste the file, or open the `.py` file from the script editor.
3. Save the script inside Mecway if you want it in the script list.

You can also clone this repository and point Mecway at `AutoBondedContacts.py`.

## Run

1. Open the model and save it (**File → Save**). The script refuses to continue if the model has never been saved.
2. Run the script.
3. Read the usage dialog and confirm it.
4. If the mesh is not usable, a second dialog explains why and the script stops. No tolerance is asked.
5. Enter the position tolerance in millimetres, for example `0.5`. Comma or dot are both accepted.
6. Enter the maximum inclination in degrees. Press OK on an empty field to keep the default of **15**. Allowed range is 0 to 90.
7. Leave the progress window in front. It has no close box and closes itself when the script finishes.
8. Read the summary dialog.
9. **Reload the model. Do not save before the reload.**

On disk the script writes:

- `yourmodel.liml.bak` — the file as Mecway saved it, before contacts were inserted.
- `yourmodel.liml` — the same file with the new face selections and bonded contacts.

If a previous run created selections named `AutoM_*`, `AutoS_*` or `AutoFace_*`, those selections and the contacts that reference them are removed and written again.

## Dialogs

All dialogs are English.

| Step | What you see |
| --- | --- |
| Start | Short description of the script, the two inputs, the node-projection rule, and the reload warning. |
| Bad mesh | `No suitable mesh.` plus the reason: no elements, only beams, or fewer than two solid/shell components. |
| Tolerance | `Position tolerance in mm (e.g. 0.5):` |
| Angle | `Max. face inclination in degrees [15]:` |
| Progress | Status text and a bar: reading components, free faces, contact pairs, coverage, node assignment, save, write. |
| Result | How many contacts were written, which component is master and slave, how many faces were dropped, and the backup path. |

Cancel either input and the script stops without changing the file. The mesh check happens before the inputs, so a bad mesh never asks for numbers.

## How contacts are chosen

Every unordered pair of components is considered. Components whose bounding boxes, expanded by the tolerance, do not meet are skipped.

### 1. Free faces

For each solid element the script asks Mecway for every face (`mw.face_nodes`). A face is free when no other element in the same component has the same set of node ids. Shared internal faces appear twice and are discarded. Shells contribute one surface per element (face 5, or face 6 if face 5 has no nodes).

### 2. Seed pairs

Two free faces, one on each component, become a seed when all of the following are true:

- The angle between their planes is at most the value you entered. The test is `|n1 · n2| >= cos(angle)`. Opposite normals still match, because only the absolute value is used. 15° is about `|n1 · n2| >= 0.9659`.
- Their bounding boxes overlap within the tolerance.
- One of its nodes is within the tolerance of an opposite face, measured along that face normal. Every free face that uses this node is then checked, not only the first. The node does not have to land inside the face outline.
- The two faces have a projected area overlap. Every opposite face that overlaps the free face of such a node is kept. A shared edge has no area and is not enough.

### 3. Coverage

After the seeds, every remaining free face is tested the same way: a node within the normal tolerance, then a projected overlap with a partner face.

### 4. One node, one contact

A node must not appear in more than one contact set. The same component can mate with several partners, and those interfaces share edge nodes.

Each accepted face is scored against its partner:

- primary: fraction of the face nodes that project onto the partner;
- secondary: smaller gap between the two planes wins ties.

Faces are assigned in score order. If any node of a face is already owned by another contact, the whole face is dropped from this contact. It cannot be split. If either side of a contact becomes empty, that contact is discarded and reported in the summary.

The better geometric fit keeps the face. A large contact does not automatically steal the last pair.

### 5. Master and slave

The side with more faces is the slave (`AutoS_n`). The other side is the master (`AutoM_n`). If both sides have the same count, the first component of the pair is the master. Mecway projects the slave onto the master, so the finer side should be the slave.

Shell face ids are chosen after pairing: face 5 if that normal points toward the partner centroid, otherwise face 6. Both sides use the same shell nodes.

## What is written

Selections and contacts are inserted immediately before the `<solution>` element.

```xml
<faceselection name="AutoM_1">
  <face eid="1" faceid="4" />
</faceselection>
<faceselection name="AutoS_1">
  <face eid="2" faceid="6" />
</faceselection>
<contact2 selectionmaster="AutoM_1" selectionslave="AutoS_1" type="bondedce" selfcontact="false">
  <tolerance>0.5 mm</tolerance>
  <stiffnessfactor>1</stiffnessfactor>
  <normalstiffness>0 GPa/m</normalstiffness>
  <tangentialstiffness>0 GPa/m</tangentialstiffness>
  <frictioncoefficient>0</frictioncoefficient>
  <conductanceperarea>0 W/m²/K</conductanceperarea>
  <tolerancetype>absolute</tolerancetype>
</contact2>
```

`type="bondedce"` is Mecway's bonded contact (constraint equations). Stiffness, friction and conductance are written as zeros, matching an otherwise empty bonded pair. The position tolerance is the value you entered, with `tolerancetype` `absolute`.

Check the contact in the Mecway tree after reload. If a build of Mecway does not accept `absolute`, set the tolerance type once in the UI and keep the face selections.

## Progress window

The window is a .NET `System.Windows.Forms` dialog (`TopMost`, no close box) updated with `Application.DoEvents`, so Mecway's UI thread can paint while the search runs. Status lines name the component and the pair index. Updates are throttled to about 12 per second.

If WinForms cannot be loaded, the search still runs; you only lose the live bar.

## Failure and recovery

| Situation | Result |
| --- | --- |
| Model has no path | Dialog `Please save the model first.` The `.liml` is not edited. |
| No seed inside the tolerance and angle | Dialog, no write. The `.bak` is not created because the file is not modified. |
| Script crashes after `os.remove` of the `.liml` | Restore `yourmodel.liml.bak`. The write goes to `yourmodel.liml.autotie.tmp` first; the original is removed only after that temp file exists. |
| You save in Mecway before reloading | Mecway writes the in-memory model, which does not contain the new contacts, and the edit is lost. Reload first. |
| A node sits on two real interfaces | The face stays with the contact that covers it better. The summary counts removed faces and dropped sets. |

## Limits

- One bonded contact per component pair, not one contact per disconnected pad.
- Shell edges and beam ends are never contacts.
- Plane elements that Mecway reports as `tri3` or `quad4` are treated as shells.
- The angle is measured between each candidate face and the partner face. A face past the limit is ignored even when a neighbour is still inside it.
- Very large tolerances pull in every face whose projection still meets the partner inside that gap.
- A node within the tolerance past an edge still pulls that face into the contact. A node farther outside does not.
- A pure crossing, where the areas overlap but no node projects onto the other face, is not a contact.
- Detection is geometric. It does not read existing Mecway contacts, named selections, or material ids.
- The script is IronPython 2.7 inside Mecway and CPython 3 outside Mecway. Outside Mecway, `python AutoBondedContacts.py` runs the geometry self-test and does not need a model.

## Self-test

From a normal Python 3 interpreter:

```bash
python AutoBondedContacts.py
```

Expected output:

```text
self-test ok
```

The test checks parallel faces, a 10° pass and a 20° reject, edge-only rejection, a crossing pair whose nodes lie outside each other, a 30-face chain, and exclusive node assignment by projected-area score.

## Constants

Edit these near the top of `AutoBondedContacts.py` if the defaults are wrong for a project:

| Name | Default | Meaning |
| --- | --- | --- |
| `CONTACT_TYPE` | `bondedce` | Mecway contact type written into `<contact2>`. |
| `DEFAULT_ANGLE_DEG` | `15` | Shown in the angle prompt. Empty input keeps this value. |

## Version

1.0.7 — a bend stays one contact. Growth follows the surface and adds the facing partner face, not only faces already marked.
1.0.6 — gaps in a contact patch are filled so the marked surface stays continuous.
1.0.5 — every free face of a qualifying node is checked, not only the first.
1.0.4 — a node counts by normal distance. Opposite faces are those with a projected overlap.
1.0.3 — contact finding uses node projection. A face is added when a node lands on the partner within the tolerance, including just past the edge.
1.0.2 — broadphase grid. Each face is tested only against nearby faces, and the overlap area is computed once.
1.0.1 — coverage is the projected polygon intersection, not a node count. A face is added as soon as that area is real.
1.0.0 — first release. Solids and shells, angle limit, progress window, exclusive nodes, `.liml` bonded contacts.
