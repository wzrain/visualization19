#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkCoordinate.h>
#include <vtkSphereSource.h>
#include <vtkButtonWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>

#include <vtkSmartPointer.h>
#include <vtkMetaImageReader.h>
#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>
#include <vtkOpenGLGPUVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkCamera.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetMapper.h>

#include <iostream>
#include <string>

void readVTKFile(const char *fileName,const char *attribute)
{
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	if(std::string(fileName) == std::string("")){reader->SetFileName("../../../clams/CLaMS_Nabro.vtk");}
	else{reader->SetFileName(fileName);}
	reader->Update();

	vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
	// cout << data->GetPointData()->GetArrayName(1) << endl;
	if(std::string(attribute) == std::string("")){
		auto time = data->GetPointData()->GetArray("time");
		auto altitude = data->GetPointData()->GetArray("altitude");
		auto pressure = data->GetPointData()->GetArray("pressure");
		auto temperature = data->GetPointData()->GetArray("temperature");
		auto potTemperature = data->GetPointData()->GetArray("pot_temperature");
		auto potVorticity = data->GetPointData()->GetArray("pot_vorticity");
		auto seedId = data->GetPointData()->GetArray("seedId");
		// std::cout << "time:  type - " << typeid(time).name() << "    elementNum - " << time->GetElementComponentSize() << std::endl;
		std::cout << "altitude:  type - " << typeid(altitude).name() << "    elementNum - " << altitude->GetNumberOfTuples() << std::endl;
		std::cout << altitude->GetTuple1(0) << std::endl;
		// std::cout << "pressure:  type - " << typeid(pressure).name() << "    elementNum - " << pressure->GetElementComponentSize() << std::endl;
		// std::cout << "temperature:  type - " << typeid(temperature).name() << "    elementNum - " << temperature->GetElementComponentSize() << std::endl;
		// std::cout << "potTemperature:  type - " << typeid(potTemperature).name() << "    elementNum - " << potTemperature->GetElementComponentSize() << std::endl;
		// std::cout << "potVorticity:  type - " << typeid(potVorticity).name() << "    elementNum - " << potVorticity->GetElementComponentSize() << std::endl;
		// std::cout << "seedId:  type - " << typeid(seedId).name() << "    elementNum - " << seedId->GetElementComponentSize() << std::endl;
		// data->GetPointData()->SetScalars(time);
		data->GetPointData()->SetScalars(altitude);
		// data->GetPointData()->SetScalars(pressure);
		// data->GetPointData()->SetScalars(temperature);
		// data->GetPointData()->SetScalars(potTemperature);
		// data->GetPointData()->SetScalars(potVorticity);
		// data->GetPointData()->SetScalars(seedId);
	}
	else{
		auto ashIndex = data->GetPointData()->GetArray(attribute);
		data->GetPointData()->SetScalars(ashIndex);
	}
	

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(reader->GetOutputPort());
	//mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(0, 0, 1);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();

	renderer->AddActor(actor);
	renderer->SetActiveCamera(aCamera);
	renderer->ResetCamera();
	renderer->SetBackground(.5, .5, .5);
	renderer->ResetCameraClippingRange();

	renWin->SetSize(600, 600);
	renWin->Render();
	renWinInteractor->Initialize();
	renWinInteractor->Start();
}

int main(int argc, char *argv[]){
	if(argc == 3){readVTKFile(argv[1],argv[2]);}
	else{readVTKFile("","");}
	return 0;
}
