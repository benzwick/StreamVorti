"""ParaView macro: Velocity vector glyphs for StreamVorti output.

MFEM writes 2D velocity as a 2-component vector ("Velocity_0", "Velocity_1").
ParaView requires a 3-component vector for glyphs. This macro creates a
Calculator filter to construct the 3D vector, then applies arrow glyphs.

Usage:
  1. Open the .pvd file in ParaView
  2. Macros → Add new macro → select this file
  3. Click the macro button to apply

Works with any StreamVorti simulation output (cavity, cylinder, channel, etc.).
"""

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

# Find the active source (the .pvd file the user opened)
source = GetActiveSource()
if source is None:
    print("ERROR: No active source. Open a .pvd file first.")
    raise SystemExit

# Create Calculator to combine 2-component velocity into 3D vector
calc = Calculator(registrationName='Velocity3D', Input=source)
calc.ResultArrayName = 'Velocity3D'
calc.Function = '"Velocity_0"*iHat + "Velocity_1"*jHat'

# Show in active view
view = GetActiveViewOrCreate('RenderView')
calcDisplay = Show(calc, view, 'UnstructuredGridRepresentation')
calcDisplay.Representation = 'Surface'
Hide(source, view)

# Color by velocity magnitude
ColorBy(calcDisplay, ('POINTS', 'Velocity3D', 'Magnitude'))
calcDisplay.RescaleTransferFunctionToDataRange(True, False)
calcDisplay.SetScalarBarVisibility(view, True)

# Create arrow glyphs
glyph = Glyph(registrationName='VelocityGlyphs', Input=calc,
              GlyphType='Arrow')
glyph.OrientationArray = ['POINTS', 'Velocity3D']
glyph.ScaleArray = ['POINTS', 'Velocity3D']
glyph.GlyphMode = 'All Points'
glyph.ScaleFactor = 0.25

# Show glyphs
glyphDisplay = Show(glyph, view, 'GeometryRepresentation')
glyphDisplay.Representation = 'Surface'
glyphDisplay.SetScalarBarVisibility(view, True)

# Set 2D camera
view.InteractionMode = '2D'
view.Update()

RenderAllViews()
