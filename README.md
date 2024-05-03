To run a simulation, create grid with StGridV2, add obstacles to the grid, generate its operators with grid.GenerateOperators (first obstacles, then operators, the order is important) then use the simulation from the Run class to run a simulation. 

The output is a VTK file in the data folder readable on Paraview.

/!\ StGrid and FlowCpp are obsolet and should not be used. Use StGridV2 and FlowCppV2 instead.

This ReadMe is a work in progess.
