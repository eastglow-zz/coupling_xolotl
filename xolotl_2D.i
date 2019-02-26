[Mesh]
  type = PETScDMDAMesh
  XolotlInput_path_name = './params_NE_2D.txt'
[]

[AuxVariables]
  [./Auxv]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Problem]
 type = XolotlProblem
 sync_variable = Auxv
[]

[Executioner]
 type = Transient
[]
