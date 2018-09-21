set(PROJECT_DIR ${PROJECT_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/.." PARENT_SCOPE)

set(JMM_CPPLibs_Dir "${CMAKE_CURRENT_SOURCE_DIR}/../JMM_CPPLibs")

set(LinearAlgebra_Dir "${JMM_CPPLibs_Dir}/LinearAlgebra")
set(Output_Dir "${JMM_CPPLibs_Dir}/Output")
set(DataStructures_Dir "${JMM_CPPLibs_Dir}/DataStructures")
set(Macros_Dir "${JMM_CPPLibs_Dir}/Macros")

set(LinearAlgebra_Headers
	${LinearAlgebra_Dir}/AffineTransformType.h
	${LinearAlgebra_Dir}/ArrayType.h
	${LinearAlgebra_Dir}/AsymmetricQuadraticNorm.h
	${LinearAlgebra_Dir}/BasisReduction.h
	${LinearAlgebra_Dir}/BasisReduction.hpp
	${LinearAlgebra_Dir}/DifferentiationType.h
	${LinearAlgebra_Dir}/FriendOperators.h
	${LinearAlgebra_Dir}/MatrixType.h
	${LinearAlgebra_Dir}/PointBaseType.h
	${LinearAlgebra_Dir}/PointType.h
	${LinearAlgebra_Dir}/RanderNorm.h
	${LinearAlgebra_Dir}/SquareCube.h
	${LinearAlgebra_Dir}/SymmetricMatrixPair.h
	${LinearAlgebra_Dir}/SymmetricMatrixType.h
	${LinearAlgebra_Dir}/VectorPairType.h
	${LinearAlgebra_Dir}/VectorType.h
	${LinearAlgebra_Dir}/VoronoiReduction.h
	PARENT_SCOPE
	)

set(Output_Headers
	${Output_Dir}/BaseIO.h
	${Output_Dir}/EnumToString.h
	${Output_Dir}/ExceptionMacro.h
	${Output_Dir}/ExportMacros.h
	${Output_Dir}/FileIO.h
	${Output_Dir}/IO.h
	${Output_Dir}/MathematicaIO.h
	${Output_Dir}/MathematicaIOExternC.h
	${Output_Dir}/MexIO.h
	${Output_Dir}/PythonIO.h
	PARENT_SCOPE
	)

set(DataStructures_Headers
	${DataStructures_Dir}/CappedVector.h	${DataStructures_Dir}/RangeAccessor.h	${DataStructures_Dir}/RedeclareTypesMacro.h	${DataStructures_Dir}/ShallowMap.h	${DataStructures_Dir}/VirtualConstIterator.h
	PARENT_SCOPE
	)

set(Macros_Headers
	${Macros_Dir}/Expand.h	${Macros_Dir}/PPCat.h	${Macros_Dir}/String.h
	PARENT_SCOPE
	)