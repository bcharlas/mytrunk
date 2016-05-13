 
/*************************************************************************
*  Copyright (C) 2016 by Benoit Charlas <benoit.charlas@gmail.com>       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//REMESHING not working!!! --> getInfo to be defined


#ifdef YADE_CGAL
#ifdef FLOW_ENGINE

#include "FlowEngine_ThermalFlowEngineT.hpp"
#include "FlowEngine.hpp"

#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>
#include<core/Scene.hpp>

//#include<core/BodyContainer.hpp>
#include<core/Body.hpp>


#include<algorithm>

const double g_Pi = 3.14159265358979323846;

/// We can add data to the Info types by inheritance
class ThermalCellInfo : public FlowCellInfo_ThermalFlowEngineT
{
	public:
	//
	
	//private:
	double cellTemperature;
	double incTemperature;// increment of temperature
	
	double cellHeatCond;
	double cellCv;
	double cellCp;
	
	double cellMolarMass;
	
	//
	// Initializing
	ThermalCellInfo (void) : FlowCellInfo_ThermalFlowEngineT() {cellTemperature=273.;}
	inline double& temperature (void) {return cellTemperature;}
	inline const double& shiftedTemperature (void) const {return cellTemperature;} 
	inline void getInfo (const ThermalCellInfo& otherCellInfo) {FlowCellInfo_ThermalFlowEngineT::getInfo(otherCellInfo); temperature()=otherCellInfo.shiftedTemperature();}
	
};

class ThermalVertexInfo : public FlowVertexInfo_ThermalFlowEngineT {
	public:
	
	double bodyTemperature;
	//bodyTemperature=0.;
	double incTemperature;// increment of temperature
	//incTemperature=0.;
	double bodyCp;
	double bodyHeatCond;// heat conductivity of the body
	double bodyMolarMass;
	//
	bool isThermal;
	//
	bool includeThermalExpansion;
	double thermalExpCoefficient;
	//
	//
	ThermalVertexInfo (void) : FlowVertexInfo_ThermalFlowEngineT() {
	  bodyTemperature=273;
	  incTemperature=0.;
	  bodyCp=1.;
	  bodyHeatCond=1.;
	  //
	  isThermal=false;
	  includeThermalExpansion=false;
	}
	//includeThermalExpansion = false;
	//
	//void getBodyVolume();
	//void setBodyVolume();
	
	//private:
	//double bodyVolume;
	
};

typedef TemplateFlowEngine_ThermalFlowEngineT<ThermalCellInfo,ThermalVertexInfo> ThermalFlowEngineT;
REGISTER_SERIALIZABLE(ThermalFlowEngineT);
//CREATE_LOGGER(ThermalFlowEngineT );
YADE_PLUGIN((ThermalFlowEngineT));


class ThermalFlowEngine : public ThermalFlowEngineT
{
	public :
		//void initCellsProperties(double temperature, double Cv, double Cp);
		//CELL_SCALAR_SETTER(double, .cellTemperature, setCellTemperature);
		void setCellTemperature(int id, double temperature);
		double getCellTemperature(int id);
		//
		void setCellCv(int id, double temperature);
		Real getCellCv(int id);
		//
		void setCellCp(int id, double temperature);
		double getCellCp(int id);
		//
		void initCellsProperties(double temperature, double cv, double cp);
		//
		//
		void setBodyTemperature(int id, double temperature);
		double getBodyTemperature(int id);
		//
		void setBodyCp(int id, double cp);
		double getBodyCp(int id);
		//
		void initBodiesProperties(boost::python::list listBodies, double temperature, double cp, double heatCond);
		//
		//void initInc(double bodyInc, double cellInc);
		//
		//void thermalExchange(int exchange);
		//void thermalExchange();
		//
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalFlowEngine,ThermalFlowEngineT,"This engine manages thermal exchanges between solids as well as between solid and gas.",
		((double, bodyIniTemperature, 273,,"Initial temperature of the bodies in Kelvin"))
		((double, bodyMolarMass, 1,,"Molar mass of the bodies in g/ mol MANQUE Unit"))
		((double, bodyCp, 1,,"Cp of the bodies in MANQUE Unit"))
		((double, bodyHeatCond, 1,,"Heat conductivity of the bodies in MANQUE Unit"))
		//((double, bodyCv, 1,,"Cv of the bodies in MANQUE Unit"))
		//
		((double, fluidMolarMass, 1,,"Molar mass of the fluid in g/ mol MANQUE Unit"))
		((double, fluidCp, 1,,"Cp of the fluid in MANQUE Unit"))
		((double, fluidCv, 1,,"Cv of the fluid in MANQUE Unit"))
		((double, fluidIniTemperature, 273,,"Initial temperature of the fluid in Kelvin"))
		((double, fluidHeatCond, 1,,"Heat conductivity of the fluid in MANQUE Unit"))
		//
		((bool,isCompressible,true,,"is the fluid compressible"))
		//
		((double, gasR, 8.314,,"Ideal gas constant, usually in J/K/mol"))
		,,,
		// cell properties
		
		.def("setCellTemperature",&ThermalFlowEngine::setCellTemperature,(boost::python::arg("id"),
										  boost::python::arg("temperature")),"set temperature properties in cell 'id'.")
		.def("getCellTemperature",&ThermalFlowEngine::getCellTemperature,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
		//
		
		.def("setCellCv",&ThermalFlowEngine::setCellCv,(boost::python::arg("id"),
								boost::python::arg("cellCv")),"set thermal properties in cell 'id'.")
		.def("getCellCv",&ThermalFlowEngine::getCellCv,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
		//
		.def("setCellCp",&ThermalFlowEngine::setCellCp,(boost::python::arg("id"),
								boost::python::arg("cellCp")),"set thermal properties in cell 'id'.")
		.def("getCellCp",&ThermalFlowEngine::getCellCp,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
		//
		.def("initCellsProperties",&ThermalFlowEngine::initCellsProperties,(boost::python::arg("temperature"),
										    boost::python::arg("cv"), 
										    boost::python::arg("cp")),"set thermal properties in cells 'id'.")
		//
		//
		
		// Bodies properties
		.def("setBodyTemperature",&ThermalFlowEngine::setBodyTemperature,(boost::python::arg("id"),
										  boost::python::arg("bodyTemperature")),"set temperature properties of the body 'id'.")
		.def("getBodyTemperature",&ThermalFlowEngine::getBodyTemperature,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
		//
		.def("setBodyCp",&ThermalFlowEngine::setBodyCp,(boost::python::arg("id"),boost::python::arg("bodyCp")),"set temperature of body 'id'.")
		.def("getBodyCp",&ThermalFlowEngine::getBodyCp,(boost::python::arg("id")),"get thermal properties of body 'id'.")
		//
		//
		
		.def("initBodiesProperties",&ThermalFlowEngine::initBodiesProperties,(boost::python::arg("listBodies"), 
										      boost::python::arg("temperature"),
										      boost::python::arg("cp"),
										      boost::python::arg("heatCond"))
										      ,"set initial thermal properties in vertices.")
		//.def("thermalExchange",&ThermalFlowEngine::thermalExchange,(boost::python::arg("exchange")),"actualize the heat exchange in the model")
		/*.def("initInc",&ThermalFlowEngine::initInc,(boost::python::arg("bodyInc"),
							    boost::python::arg("cellInc")),
							    "initialize the temperature increments")
		*/
		//.def("setBodyVolume",&ThermalFlowEngine::setCellTemperature,(boost::python::arg("id"),boost::python::arg("bodyVolume")),"set volume of body 'id'.")
		//.def("getBodyTemperature",&ThermalFlowEngine::getCellTemperature,(boost::python::arg("id")),"get temperature in cell 'id'.")
		)
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ThermalFlowEngine);



void ThermalFlowEngine::setCellTemperature(int id, double temperature)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellTemperature = temperature;
}

double ThermalFlowEngine::getCellTemperature(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellTemperature;
}
//

void ThermalFlowEngine::setCellCv(int id, double Cv)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellCv = Cv;
}

double ThermalFlowEngine::getCellCv(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellCv;
}
//
void ThermalFlowEngine::setCellCp(int id, double Cp)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellCp = Cp;
}

double ThermalFlowEngine::getCellCp(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellCp;
}
//
void ThermalFlowEngine::initCellsProperties(double temperature, double cv, double cp)
{
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	 cell->info().cellTemperature = temperature;
	 cell->info().cellCv = cv;
	 cell->info().cellCp = cp;
	 fluidBulkModulus=   gasR/cell->info().cellMolarMass * temperature;// MANQUE mass / volume
	}
}

//


void ThermalFlowEngine::setBodyTemperature(int id, double temperature)
{
	solver->T[solver->currentTes].vertexHandles[id]->info().bodyTemperature = temperature;
}

double ThermalFlowEngine::getBodyTemperature(int id)
{
	return solver->T[solver->currentTes].vertexHandles[id]->info().bodyTemperature;
}
//
void ThermalFlowEngine::setBodyCp(int id, double cp)
{
	solver->T[solver->currentTes].vertexHandles[id]->info().bodyCp = cp;
}

double ThermalFlowEngine::getBodyCp(int id)
{
	return solver->T[solver->currentTes].vertexHandles[id]->info().bodyCp;
}
//
//
void ThermalFlowEngine::initBodiesProperties(boost::python::list listBodies, double temperature, double cp, double heatCond)
{
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	  int id=cell->info().id;
	  if( listBodies.count(id) >0 )
	  //if(found )
	  {
	    solver->T[solver->currentTes].vertexHandles[id]->info().isThermal=true;
	    solver->T[solver->currentTes].vertexHandles[id]->info().bodyTemperature=temperature;
	    solver->T[solver->currentTes].vertexHandles[id]->info().bodyCp=cp;
	    solver->T[solver->currentTes].vertexHandles[id]->info().bodyHeatCond=heatCond;
	  }
	}
}
//solver->T[solver->currentTes].cellHandles[id]->vertex(i)->info().id()

/*
void ThermalFlowEngine::initInc(double bodyInc, double cellInc)
{
  //Initializing Body temperature increment
  FOREACH(VertexHandle& vertex, solver->T[solver->currentTes].vertexHandles)
  {
    vertex->info().incTemperature=bodyInc;
  }
  //Initializing Cell temperature increment
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
  {
    cell->info().incTemperature=cellInc;
  }
}
*/

/*
void ThermalFlowEngine::thermalExchange(int exchange)
{
  //shared_ptr<BodyContainer> bodies;
  
  //Initializing Body temperature increment
  FOREACH(VertexHandle& vertex, solver->T[solver->currentTes].vertexHandles)
  {
    vertex->info().incTemperature=0.;
  }
  //Initializing Cell temperature increment
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
  {
    cell->info().incTemperature=0.;
  }
  
  // heat conductivity between solids
  FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
  {
    if(!I->isReal()) continue;
    ScGeom*    geom= static_cast<ScGeom*>(I->geom.get());
    
    //if(!geom->radius1 || !geom->radius2) continue;// if one of the bodies is not a sphere - potentially to remove
    
    if( !solver->T[solver->currentTes].vertexHandles[I->getId1()]->info().isThermal 
      || !solver->T[solver->currentTes].vertexHandles[I->getId2()]->info().isThermal) continue;
    
    double R1=geom->radius1;
    double R2=geom->radius2;
    
    double L1=solver->T[solver->currentTes].vertexHandles[I->getId1()]->info().bodyHeatCond;
    double L2=solver->T[solver->currentTes].vertexHandles[I->getId2()]->info().bodyHeatCond;
    
    double T1=solver->T[solver->currentTes].vertexHandles[I->getId1()]->info().bodyTemperature;
    double T2=solver->T[solver->currentTes].vertexHandles[I->getId2()]->info().bodyTemperature;
    
    
    double& un=geom->penetrationDepth;
    
    double Rmin=std::min(R1,R2);
    
    solver->T[solver->currentTes].vertexHandles[I->getId1()]->info().incTemperature +=-L1*L2/((R1-un/2.)*L2 + (R2-un/2.)*L1) * g_Pi *Rmin * un*(T1-T2) * scene->dt;
    solver->T[solver->currentTes].vertexHandles[I->getId2()]->info().incTemperature += L1*L2/((R1-un/2.)*L2 + (R2-un/2.)*L1) * g_Pi *Rmin * un *(T1-T2) *scene->dt;
    
    
  }
  
  // heat exchange with neighboring bodies and fluid cells
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
  {
    for (unsigned int ngb=0;ngb<4;ngb++)
    {
      double Ti = cell->info().cellTemperature;
      double Tb = cell->vertex(ngb)->info().bodyTemperature;
      double Tc = cell->neighbor(ngb)->info().cellTemperature;
      
      int bodyId = cell->vertex(ngb)->info().id();
      // WARNING doesn't work if the considered bodies are not spheres
      const Sphere* sphere = YADE_CAST<Sphere*> (Body::byId(bodyId)->shape.get());
      // manque rayon
      double R=sphere->radius;
      
      // density 
      double rho=cell->info().p() * cell->info().cellMolarMass / (cell->info().cellTemperature * gasR);
      
      // solid / gas - WARNING approximative estimation of the surface
      cell->info().incTemperature = -(1.14+0.85/R)*g_Pi*pow(R,2)*pow((Ti - Tb),1.25) *scene->dt / cell->info().cellCp /rho;
      cell->vertex(ngb)->info().incTemperature = (1.14+0.85/R)*g_Pi*pow(R,2)*pow((Ti - Tb),1.25) *scene->dt / cell->vertex(ngb)->info().bodyCp * 4/3 * g_Pi * pow(R,3) / Body::byId(bodyId)->state->mass ;
      
      //gas / gas - all the gas is supposed to have the same heat conductivity
      cell->info().incTemperature = - cell->info().cellHeatCond * (Ti - Tc);
      cell->neighbor(ngb)->info().incTemperature = - cell->info().cellHeatCond * (Ti - Tc);
      
    }
  }
  
  // incrementing the temperature
  FOREACH(VertexHandle& vertex, solver->T[solver->currentTes].vertexHandles)
  {
    vertex->info().bodyTemperature +=vertex->info().incTemperature;
    //
    vertex->info().incTemperature = 0;
  }
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
  {
    // density of the gas in the cell
    //double rho = cell->info().p() * cell->info().cellMolarMass / gasR / cell->info().cellTemperature
    //
    
    if(isCompressible)
      {
	// pressure change due to temperature
	cell->info().p() -=cell->info().p()  * ( cell->info().incTemperature / cell->info().cellTemperature ) *scene->dt;
	// pressure change due to mass flow between cells
	for (unsigned int ngb=0;ngb<4;ngb++)
	{
	  cell->info().p() -= gasR / cell->info().cellMolarMass / getCellVoidVolume(cell->info().id) * ( cell->info().kNorm())[ngb]*(cell->info().p()-cell->neighbor(ngb)->info().p()) *scene->dt;
	}
	
	
      }
    cell->info().cellTemperature +=cell->info().incTemperature;
    
  }
  
}

*/


//YADE_PLUGIN((ThermalFlowEngineT));
CREATE_LOGGER(ThermalFlowEngine );
YADE_PLUGIN((ThermalFlowEngine));

#endif //FLOW_ENGINE

#endif /* YADE_CGAL */