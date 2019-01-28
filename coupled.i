[Mesh]
 type = GeneratedMesh
 dim = 2
 nx = 20
 ny = 10
[]

[Variables]
    [./u]
    [../]
[]

[AuxVariables]
    [./v]
    [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./td]
    type = TimeDerivative
    variable = u
  [../]
  [./cf]
    type = CoupledForce
    coef = 1
    variable = u
    v=v
  [../]
[]

[BCs]
    [./left]
        type = DirichletBC
        variable = u
        boundary = left
        value = 0
    [../]
    [./right]
        type = DirichletBC
        variable = u
        boundary = right
        value = 0
    [../]
[]

[Executioner]
 type = Transient
 num_steps = 10
 dt = 10000

 solve_type = 'PJFNK'

 petsc_options_iname = '-pc_type -pc_hypre_type'
 petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
 exodus = true
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        input_files = 'xolotl.i'
        app_type = coupling_xolotlApp
        library_path = 'lib'
    [../]
[]

[Transfers]
    [./fromsub]
        type = MultiAppNearestNodeTransfer
        direction = from_multiapp
        multi_app = sub_app
        source_variable = u
        variable = v
        fixed_meshes = true
    [../]
[]
