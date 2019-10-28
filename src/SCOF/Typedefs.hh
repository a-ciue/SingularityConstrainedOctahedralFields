/*
 * Typedefs.hh
 *
 *  Created on: May 11, 2018
 *      Author: hliu
 */
#pragma once

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


namespace OVM = OpenVolumeMesh;
namespace OM = OpenMesh;

namespace SCOF
{
//OVM
	using VH = OVM::VertexHandle;
	using HEH = OVM::HalfEdgeHandle;
	using EH = OVM::EdgeHandle;
	using HFH = OVM::HalfFaceHandle;
	using FH = OVM::FaceHandle;
	using CH = OVM::CellHandle;

	template <typename T>
	using VP = OVM::VertexPropertyT<T>;
	template <typename T>
	using HEP = OVM::HalfEdgePropertyT<T>;
	template <typename T>
	using EP = OVM::EdgePropertyT<T>;
	template <typename T>
	using HFP = OVM::HalfFacePropertyT<T>;
	template <typename T>
	using FP = OVM::FacePropertyT<T>;
	template <typename T>
	using CP = OVM::CellPropertyT<T>;
	template <typename T>
	using MP = OVM::MeshPropertyT<T>;

//OM
	struct TriTraits : public OpenMesh::DefaultTraits
	{
	  typedef OVM::Geometry::Vec3d Point;
	};
	typedef OpenMesh::TriMesh_ArrayKernelT<TriTraits>  TriMesh;

	using TriVH = OM::VertexHandle;
	using TriHEH = OM::HalfedgeHandle;
	using TriEH = OM::EdgeHandle;
	using TriFH = OM::FaceHandle;

	template <typename T>
	using TriVP = OM::VPropHandleT<T>;
	template <typename T>
	using TriHP = OM::HPropHandleT<T>;
	template <typename T>
	using TriEP = OM::EPropHandleT<T>;
	template <typename T>
	using TriFP = OM::FPropHandleT<T>;
}
