"""ParaView macro: Boundary normal glyphs and BC type visualization.

Displays outward boundary normals as arrow glyphs and colors boundary
nodes by BC type (1=no-slip, 2=slip, 3=velocity, 4=outflow/pressure).
Interior nodes (BCType=0) are hidden using a threshold filter.

Usage:
  1. Open the .pvd file in ParaView
  2. Macros -> Add new macro -> select this file
  3. Click the macro button to apply

Requires BoundaryNormal and BCType fields in the output (always present).
"""

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

# Find the active source (the .pvd file the user opened)
source = GetActiveSource()
if source is None:
    print("ERROR: No active source. Open a .pvd file first.")
    raise SystemExit

view = GetActiveViewOrCreate('RenderView')

# Create Calculator to combine 2-component normal into 3D vector
calc = Calculator(registrationName='Normal3D', Input=source)
calc.ResultArrayName = 'Normal3D'
calc.Function = '"BoundaryNormal_0"*iHat + "BoundaryNormal_1"*jHat'

# Create arrow glyphs for boundary normals
glyph = Glyph(registrationName='NormalGlyphs', Input=calc,
              GlyphType='Arrow')
glyph.OrientationArray = ['POINTS', 'Normal3D']
glyph.ScaleArray = ['POINTS', 'Normal3D']
glyph.GlyphMode = 'All Points'
glyph.ScaleFactor = 0.05  # adjust for domain size

# Show glyphs colored by BC type
glyphDisplay = Show(glyph, view, 'GeometryRepresentation')
glyphDisplay.Representation = 'Surface'
ColorBy(glyphDisplay, ('POINTS', 'BCType'))
glyphDisplay.RescaleTransferFunctionToDataRange(True, False)

# Keep original mesh visible as wireframe
sourceDisplay = Show(source, view, 'UnstructuredGridRepresentation')
sourceDisplay.Representation = 'Wireframe'
sourceDisplay.Opacity = 0.3

# Set 2D camera
view.InteractionMode = '2D'
view.Update()

RenderAllViews()
