[Mesh]
  type = PETScDMDAMesh
  XolotlInput_path_name = './params_NE_2D.txt'
[]

[AuxVariables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Problem]
 type = XolotlProblem
 sync_variable = u
[]

[Executioner]
 type = Transient
[]
