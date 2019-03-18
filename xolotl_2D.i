[Mesh]
  type = PETScDMDAMesh
  dim = 2
  XolotlInput_path_name = './params_NE_2D.txt'
[]

[AuxVariables]
  [./Auxv]
    order = FIRST
    family = LAGRANGE
    outputs = exodus
  [../]
  [./AuxGB]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Problem]
 type = XolotlProblem
 sync_variable = Auxv
 sync_GB = AuxGB
 outputs = exodus
[]

[Functions]
  [./dts]
    type = PiecewiseLinear
    x = '0 1 10 100 1000 10000 100000000'
    y = '0.1 0.1 1 10 100 1000 5000'
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    # interval = 10
    interval = 200
  [../]
[]

[Executioner]
 type = Transient
  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 0.01
  [../]
[]
