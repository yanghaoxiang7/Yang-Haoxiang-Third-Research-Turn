#ifndef __MultiCopterDriver_h__
#define __MultiCopterDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLParticles.h"
#include "Particles.h"
#include "MultiCopter.h"

template<int d> class MultiCopterDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	MultiCopter<d> copter;
	Array<VectorD> targets = {
		VectorD(0.0, 0.0, -1.0),
		VectorD(1.0, 1.0, -1.2),
		VectorD(2.0, 0.0, -1.4),
		VectorD(1.0, -1.0, -1.2),
		VectorD(0.0, 0.0, -1.0),
		VectorD(0.0, 0.0, 0.0)};

	OpenGLSegmentMesh* opengl_copter=nullptr;
	OpenGLSegmentMesh* opengl_circles=nullptr;
	Array<VectorD> circle_vtx;
	OpenGLParticles<Particles<3> >* opengl_targets=nullptr;

public:
	virtual void Initialize(const int flag)
	{
		copter.Initialize(flag);

		const VectorD dir = -VectorD::UnitZ();
		const real d1 = copter.arm_length / (real)(std::sqrt(2.0));
		copter.Add_Rotor(VectorD(d1, -d1, 0.0), dir);
		copter.Add_Rotor(VectorD(d1, d1, 0.0), dir);
		copter.Add_Rotor(VectorD(-d1, d1, 0.0), dir);
		copter.Add_Rotor(VectorD(-d1, -d1, 0.0), dir);

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		opengl_copter=Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_copter->mesh.elements.resize(2);
		opengl_copter->mesh.elements[0]=Vector2i(0,2);
		opengl_copter->mesh.elements[1]=Vector2i(1,3);
		opengl_copter->line_width=2.f;
		*opengl_copter->mesh.vertices=copter.body_rotor_pos;

		opengl_copter->Set_Data_Refreshed();
		opengl_copter->Initialize();

		opengl_circles=Add_Interactive_Object<OpenGLSegmentMesh>();
		SegmentMesh<3>& mesh=opengl_circles->mesh;
		
		for(int i=0;i<copter.body_rotor_pos.size();i++){
			Vector3 center=copter.body_rotor_pos[i];
			real r=.02f;int n=16;
			real theta=3.1415927f*2.f/(real)n;
			int start=(int)(*mesh.vertices).size();
			for(int j=0;j<n;j++){
				real angle=(real)j*theta;
				Vector3 p=center+Vector3(r*cos(angle),r*sin(angle),-0.005);
				circle_vtx.push_back(p);
				(*mesh.vertices).push_back(copter.World_Coord(p));}
			for(int j=0;j<n-1;j++){
				mesh.elements.push_back(Vector2i(start+j,start+j+1));}
			mesh.elements.push_back(Vector2i(start+n-1,start));}

		opengl_circles->Set_Color(OpenGLColor(0.f,1.f,0.f,1.f));

		opengl_circles->Set_Data_Refreshed();
		opengl_circles->Initialize();
		
		opengl_targets=Add_Interactive_Object<OpenGLParticles<Particles<3> > >();
		opengl_targets->particles.Resize((const int)targets.size());
		for(int i=0;i<targets.size();i++){
			opengl_targets->particles.X(i)=targets[i];}
		opengl_targets->Set_Data_Refreshed();
		opengl_targets->Initialize();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		for(int i=0;i<copter.body_rotor_pos.size();i++){
			(*opengl_copter->mesh.vertices)[i]=copter.World_Coord(copter.body_rotor_pos[i]);}
		opengl_copter->Set_Data_Refreshed();

		for(int i=0;i<(*opengl_circles->mesh.vertices).size();i++){
			(*opengl_circles->mesh.vertices)[i]=copter.World_Coord(circle_vtx[i]);}
		opengl_circles->Set_Data_Refreshed();
		
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		static real t = 0.0;
		int idx = (int)(t / 5.0);
		if (idx > 5) idx = 5;
		copter.Advance(dt, targets[idx]);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
		t += dt;
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}
};
#endif
