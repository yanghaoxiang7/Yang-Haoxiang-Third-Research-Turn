//#####################################################################
// Mass spring driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __MassSpringInteractiveDriver_h__
#define __MassSpringInteractiveDriver_h__
#include <memory>
#include "Common.h"
#include "Mesh.h"
#include "Driver.h"
#include "SoftBodyMassSpring.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"

template<int d> class MassSpringInteractivDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
public:
	SoftBodyMassSpring<d> soft_body;
	const real dt=(real).02;

	////mesh data for visualization only
	std::shared_ptr<SegmentMesh<d> > vis_segment_mesh=nullptr;
	std::shared_ptr<TriangleMesh<d> > vis_triangle_mesh=nullptr;

	////visualization data
	OpenGLSegmentMesh* opengl_segments=nullptr;							////vector field
	Array<OpenGLSphere*> opengl_spheres;								////spheres

	virtual void Initialize(){OpenGLViewer::Initialize();}
	virtual void Run(){OpenGLViewer::Run();}
	virtual void Initialize_Data()
	{
		Initialize_Simulation_Data();
		Initialize_OpenGL_Data();
	}

	void Initialize_OpenGL_Data()
	{		
		////initialize a segment mesh to visualize the trace
		opengl_segments=Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_segments->mesh.Vertices().resize(soft_body.particles.Size());
		for(int i=0;i<soft_body.particles.Size();i++){
			opengl_segments->mesh.Vertices()[i]=soft_body.particles.X(i);}
		opengl_segments->mesh.Elements().resize(soft_body.springs.size());
		for(int i=0;i<soft_body.springs.size();i++){
			opengl_segments->mesh.Elements()[i]=soft_body.springs[i];}
		opengl_segments->Set_Data_Refreshed();
		opengl_segments->Initialize();	

		////initialize a sphere to visualize the particle
		for(int i=0;i<soft_body.particles.Size();i++){
			OpenGLSphere* opengl_sphere=Add_Interactive_Object<OpenGLSphere>();
			opengl_sphere->pos=soft_body.particles.X(i);
			opengl_sphere->radius=(real).02;
			Set_Color(opengl_sphere,OpenGLColor(.0,1.,.0,1.));
			opengl_sphere->Set_Data_Refreshed();
			opengl_sphere->polygon_mode=PolygonMode::Fill;
			opengl_sphere->Initialize();
			opengl_spheres.push_back(opengl_sphere);}

		////set OpenGL rendering environments
		auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
		OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
		OpenGLUbos::Update_Lights_Ubo();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		////update and sync data for segments
		for(int i=0;i<opengl_segments->mesh.Vertices().size();i++){
			opengl_segments->mesh.Vertices()[i]=soft_body.particles.X(i);}
		opengl_segments->Set_Data_Refreshed();

		////update and sync data for spheres
		for(int i=0;i<opengl_spheres.size();i++){
			opengl_spheres[i]->pos=soft_body.particles.X(i);
			opengl_spheres[i]->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		Advance_One_Time_Step(dt,time+dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		soft_body.Advance(dt);
	}
	
	virtual void Initialize_Simulation_Data()
	{
		switch(test){
		case 1:{	////rod, for both 2D and 3D
			////initialize spring vertices
			real length=(real)1;int n=8;real dx=length/(real)n;
			soft_body.particles.Resize(n);
			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=VectorD::Unit(0)*(real)i*dx;
				soft_body.particles.M(i)=(real)1;}
			////initialize springs
			for(int i=0;i<n-1;i++){Vector2i s(i,i+1);
				soft_body.springs.push_back(s);}
			////set boundary conditions
			soft_body.Set_Boundary_Node(0);
		}break;
		case 2:{	////cloth, for 3D only
			////create a cloth mesh
			real length=(real)1;int width=4*scale;int height=6*scale;real step=length/(real)width;
			TriangleMesh<d> cloth_mesh;
			Build_Cloth_Mesh(width,height,step,&cloth_mesh,0,2);
			int n=(int)cloth_mesh.Vertices().size();
			Array<Vector2i> edges;Get_Mesh_Edges(cloth_mesh,edges);
			
			////copy cloth mesh vertices to spring particles 
			soft_body.particles.Resize(n);
			for(int i=0;i<n;i++){
				soft_body.particles.X(i)=cloth_mesh.Vertices()[i];
				soft_body.particles.M(i)=(real)1;}
			////copy cloth mesh edges to springs
			soft_body.springs=edges;

			////set boundary conditions
			soft_body.Set_Boundary_Node(0);
			soft_body.Set_Boundary_Node(width-1);

			////set the visualization triangle mesh
			vis_triangle_mesh.reset(new TriangleMesh<d>(soft_body.particles.XPtr()));	
			vis_triangle_mesh->elements=cloth_mesh.elements;
		}break;
		case 3:{	////volumetric beam, for 3D only
			int n=4*scale;real dx=(real)1/(real)n;
			Build_Beam_Particles_And_Springs(soft_body.particles,soft_body.springs,n,dx);
			for(int i=0;i<4;i++)soft_body.Set_Boundary_Node(i);
		}break;

		//////////////////////////////////////////////////////////////////////////
		////YOUR IMPLEMENTATION (P1 TASK): create your own mass-spring simulation
		case 4:{
			/* Your implementation */
		}break;
		}

		//////set the visualization mesh
		//vis_segment_mesh.reset(new SegmentMesh<d>(soft_body.particles.XPtr()));
		//vis_segment_mesh->elements=soft_body.springs;

		soft_body.Initialize();
	}

protected:
	////Helper functions
	void Build_Cloth_Mesh(const int cell_num_0,const int cell_num_1,const real dx,TriangleMesh<3>* mesh,int axis_0=0,int axis_1=1)
	{
		mesh->elements.resize(2*(cell_num_0-1)*(cell_num_1-1));int t=0;
		for(int i=1;i<=cell_num_0-1;i++)for(int j=1;j<=cell_num_1-1;j++){ // counterclockwise node ordering
			if(i%2){mesh->elements[t++]=Vector3i(i+cell_num_0*(j-1),i+1+cell_num_0*(j-1),i+cell_num_0*j);mesh->elements[t++]=Vector3i(i+1+cell_num_0*(j-1),i+1+cell_num_0*j,i+cell_num_0*j);}
			else{mesh->elements[t++]=Vector3i(i+cell_num_0*(j-1),i+1+cell_num_0*(j-1),i+1+cell_num_0*j);mesh->elements[t++]=Vector3i(i+cell_num_0*(j-1),i+1+cell_num_0*j,i+cell_num_0*j);}}
		for(size_type i=0;i<mesh->elements.size();i++){mesh->elements[i]-=Vector3i::Ones();
		/*swap y and z*/int tmp=mesh->elements[i][1];mesh->elements[i][1]=mesh->elements[i][2];mesh->elements[i][2]=tmp;}
		for(int j=0;j<cell_num_1;j++)for(int i=0;i<cell_num_0;i++){VectorD pos=VectorD::Zero();pos[axis_0]=(real)i*dx;pos[axis_1]=(real)j*dx;mesh->Vertices().push_back(pos);}
	}

	void Get_Mesh_Edges(const TriangleMesh<3>& mesh,Array<Vector2i>& edges)
	{
		Hashset<Vector2i> edge_hashset;ArrayF<Vector2i,6> element_edges;
		for(const auto& vtx:mesh.elements){
			edge_hashset.insert(Sorted(Vector2i(vtx[0],vtx[1])));
			edge_hashset.insert(Sorted(Vector2i(vtx[1],vtx[2])));
			edge_hashset.insert(Sorted(Vector2i(vtx[2],vtx[0])));}
		for(const auto& edge:edge_hashset)edges.push_back(edge);
	}	

	void Build_Beam_Particles_And_Springs(Particles<3>& particles,Array<Vector2i>& edges,int n,real dx,Vector3 pos=Vector3::Zero())
	{
		particles.Resize(n*4);
		for(int i=0;i<particles.Size();i++){
			particles.M(i)=(real)1;}
		for(int i=0;i<n;i++){
			particles.X(i*4)=pos+Vector3(dx*(real)i,(real)0,(real)0);
			particles.X(i*4+1)=pos+Vector3(dx*(real)i,(real)0,(real)dx);
			particles.X(i*4+2)=pos+Vector3(dx*(real)i,(real)dx,(real)0);
			particles.X(i*4+3)=pos+Vector3(dx*(real)i,(real)dx,(real)dx);
			edges.push_back(Vector2i(i*4,i*4+1));
			edges.push_back(Vector2i(i*4+1,i*4+3));
			edges.push_back(Vector2i(i*4+3,i*4+2));
			edges.push_back(Vector2i(i*4+2,i*4));
			if(i<n-1){
				edges.push_back(Vector2i(i*4,i*4+4));
				edges.push_back(Vector2i(i*4+1,i*4+5));
				edges.push_back(Vector2i(i*4+2,i*4+6));
				edges.push_back(Vector2i(i*4+3,i*4+7));
				
				edges.push_back(Vector2i(i*4,i*4+7));
				edges.push_back(Vector2i(i*4+1,i*4+6));
				edges.push_back(Vector2i(i*4+2,i*4+5));
				edges.push_back(Vector2i(i*4+3,i*4+4));}}
	}

	Vector2i Sorted(const Vector2i& v){return v[0]>v[1]?v:Vector2i(v[1],v[0]);}

	////for compiling template with d=2
	void Build_Cloth_Mesh(const int cell_num_0,const int cell_num_1,const real dx,TriangleMesh<2>* mesh,int axis_0=0,int axis_1=1){}
	void Get_Mesh_Edges(const TriangleMesh<2>& mesh,Array<Vector2i>& edges){}
	void Build_Beam_Particles_And_Springs(Particles<2>& particles,Array<Vector2i>& edges,int n,real dx,Vector2 pos=Vector2::Zero()){}
};
#endif