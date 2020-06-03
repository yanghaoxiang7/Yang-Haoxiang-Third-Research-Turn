#ifndef __GridFluidDriver_h__
#define __GridFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "Particles.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"

template<int d> class GridFluidDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	GridFluid<d> fluid;
	real v_rescale=(real).05;
	bool add_particle=false;

	OpenGLSegmentMesh* opengl_vectors=nullptr;							////vector field for velocity
	OpenGLPolygon* opengl_polygon=nullptr;								////a rectangle for domain boundary
	Array<OpenGLSolidCircle*> opengl_circles;							////advective particles
	Hashset<int> invis_particles;

public:
	virtual void Initialize()
	{
		fluid.Initialize();
		fluid.Initialize_Visualization_Particles();
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		////initialize vector field (write all vectors as segments)
		opengl_vectors=Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_vectors->mesh.elements.resize(fluid.node_num);
		opengl_vectors->mesh.vertices->resize(fluid.node_num*2);

		for(int i=0;i<fluid.particles.Size();i++){
			Add_Solid_Circle(i);}

		for(int i=0;i<fluid.node_num;i++){
			opengl_vectors->mesh.elements[i]=Vector2i(i*2,i*2+1);
			(*opengl_vectors->mesh.vertices)[i*2]=V3(fluid.grid.Node(fluid.grid.Node_Coord(i)));
			VectorD pos2=fluid.grid.Node(fluid.grid.Node_Coord(i))+fluid.u[i]*v_rescale;
			(*opengl_vectors->mesh.vertices)[i*2+1]=V3(pos2);}
		opengl_vectors->Set_Data_Refreshed();
		opengl_vectors->Initialize();

		////initialize polygon
		opengl_polygon=Add_Interactive_Object<OpenGLPolygon>();
		opengl_polygon->vtx.push_back(Vector3::Zero());
		opengl_polygon->vtx.push_back(Vector3::Unit(0)*(real)2);
		opengl_polygon->vtx.push_back(Vector3::Unit(0)*(real)2+Vector3::Unit(1)*(real)1);
		opengl_polygon->vtx.push_back(Vector3::Unit(1)*(real)1);
		Set_Color(opengl_polygon,OpenGLColor(.0,1.,1.,1.));
		Set_Line_Width(opengl_polygon,4.f);
		opengl_polygon->Set_Data_Refreshed();
		opengl_polygon->Initialize();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		////sync data for vectors
		for(int i=0;i<fluid.node_num;i++){
			VectorD pos2=fluid.grid.Node(fluid.grid.Node_Coord(i))+fluid.u[i]*v_rescale;
			(*opengl_vectors->mesh.vertices)[i*2+1]=V3(pos2);}
		opengl_vectors->Set_Data_Refreshed();

		////particle visualization
		fluid.Update_Visualization_Particles(dt);
		for(int i=0;i<fluid.particles.Size();i++){
			bool outside=false;
			real epsilon=fluid.grid.dx*(real).5;
			for(int j=0;j<d;j++){
				if(fluid.particles.X(i)[j]<fluid.grid.domain_min[j]+epsilon||
					fluid.particles.X(i)[j]>fluid.grid.domain_max[j]-epsilon){outside=true;break;}}
			if(outside){
				opengl_circles[i]->visible=false;
				opengl_circles[i]->Set_Data_Refreshed();
				invis_particles.insert(i);
				continue;}

			auto opengl_circle=opengl_circles[i];
			opengl_circle->pos=V3(fluid.particles.X(i));
			opengl_circle->pos[2]=(real).05;
			opengl_circle->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////User interaction
	virtual bool Mouse_Drag(int x,int y,int w,int h)
	{
		if(!add_particle)return false;
		Vector3f win_pos=opengl_window->Project(Vector3f::Zero());
		Vector3f pos=opengl_window->Unproject(Vector3f((float)x,(float)y,win_pos[2]));
		VectorD p_pos;for(int i=0;i<d;i++)p_pos[i]=(real)pos[i];
		fluid.src_pos=p_pos;
		Add_Source_Particle(p_pos);
		return true;
	}

	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
	{
		if(left!=1&&left!=-1){return false;}
		if(left==-1&&add_particle){
			add_particle=false;
			fluid.src_pos=VectorD::Ones()*-1;	////turn off source
			return true;}

		Vector3f win_pos=opengl_window->Project(Vector3f::Zero());
		Vector3f pos=opengl_window->Unproject(Vector3f((float)x,(float)y,win_pos[2]));
		VectorD p_pos;for(int i=0;i<d;i++)p_pos[i]=(real)pos[i];
		fluid.src_pos=p_pos;
		Add_Source_Particle(p_pos);
		add_particle=true;
		return true;
	}

	////Keyboard interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('v',&Keyboard_Event_V_Func,"press v");
	}

	virtual void Keyboard_Event_V()
	{
		opengl_vectors->visible=!opengl_vectors->visible;
		opengl_vectors->Set_Data_Refreshed();
	}
	Define_Function_Object(GridFluidDriver,Keyboard_Event_V);

	void Add_Source_Particle(VectorD p_pos)
	{
		real rx=.1*static_cast<float>(rand()%1000)/1000.-.05;
		real ry=.1*static_cast<float>(rand()%1000)/1000.-.05;
		p_pos[0]+=rx;p_pos[1]+=ry;
		if(!invis_particles.empty()){
			int p=(*invis_particles.begin());
			fluid.particles.X(p)=p_pos;
			opengl_circles[p]->pos=V3(fluid.particles.X(p));
			opengl_circles[p]->Update_Model_Matrix();
			opengl_circles[p]->Set_Data_Refreshed();
			opengl_circles[p]->visible=true;
			invis_particles.erase(p);
		}
		else{
			Add_Particle(p_pos);
			Add_Solid_Circle(fluid.particles.Size()-1);	
		}
	}

	void Add_Particle(VectorD pos)
	{
		int i=fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i)=pos;
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c(1.f,0.68f,0.26f,1.f);
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circle->visible=false;
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos=V3(fluid.particles.X(i));
		opengl_circle->radius=(real).01;
		opengl_circle->color=c;
		opengl_circle->Initialize();
		opengl_circle->Update_Model_Matrix();
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->visible=true;
	}

protected:
	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif