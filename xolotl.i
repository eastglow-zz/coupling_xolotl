[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 10
 ny = 10
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
 XolotlInput_path_name = './params_NE_2D.txt'
[]

[Executioner]
 type = Transient
[]
