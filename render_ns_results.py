import sys
import paraview.simple as pvs

if __name__=="__main__":

    # Get results file
    results_file = sys.argv[-1]

    # Load into paraview
    data_reader = pvs.OpenDataFile(results_file)

    # Table to points
    points = pvs.TableToPoints(data_reader, XColumn='x', YColumn='y', ZColumn='z')
    pvs.Show(points)
    
    # Delauney 2D filter
    delaunay = pvs.Delaunay2D(points)

    # Add grid to view
    pvs.Show(delaunay)

    # Set color of grid
    point_data = points.PointData
    P_range = point_data[0].GetRange(0)
    V_range = point_data[2].GetRange(0)
    print(V_range)
    disp_dl2 = pvs.GetDisplayProperties(delaunay)
    disp_dl2.LookupTable = pvs.MakeBlueToRedLT(P_range[0], P_range[1])
    disp_dl2.ColorArrayName = 'P'

    # Initialize vector calculator
    calc = pvs.Calculator(points, Function='u*iHat+v*jHat+w*kHat')

    # Intialize arrow glyph
    arrow = pvs.Glyph(calc, ScaleArray='V', ScaleFactor=0.2/V_range[1], GlyphType='2D Glyph', OrientationArray='Result')

    # Add arrows to view
    pvs.Show(arrow)

    # Render
    pvs.Interact()
    pvs.Render()