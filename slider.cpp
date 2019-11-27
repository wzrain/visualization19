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
#include "vtkEasyTransfer.hpp"
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkProperty.h>
#include <vtkExtractVOI.h>
#include <vtkImageDataGeometryFilter.h>

#include <string>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void CreateImage(vtkSmartPointer<vtkImageData> image, unsigned char* color1, unsigned char* color2)
{
	// Specify the size of the image data
	image->SetDimensions(10, 10, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	int* dims = image->GetDimensions();

	// Fill the image with
	for (int y = 0; y < dims[1]; y++)
	{
		for (int x = 0; x < dims[0]; x++)
		{
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			if (x < 5)
			{
				pixel[0] = color1[0];
				pixel[1] = color1[1];
				pixel[2] = color1[2];
			}
			else
			{
				pixel[0] = color2[0];
				pixel[1] = color2[1];
				pixel[2] = color2[2];
			}
		}
	}
}


// Callback for the slider interaction
class vtkSliderCallback : public vtkCommand
{
public:
	static vtkSliderCallback* New() {
		return new vtkSliderCallback;
	}
	virtual void Execute(vtkObject* caller, unsigned long, void*) {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		int ind = round(static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue());
		mapper->SetInputData(data[ind]);
	}
	vtkSliderCallback() : mapper(0), data(0) {}
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vector<vtkSmartPointer<vtkPolyData>> data;
};


class vtkButtonCallback : public vtkCommand
{
public:
	static vtkButtonCallback* New() {
		return new vtkButtonCallback;
	}
	virtual void Execute(vtkObject* caller, unsigned long, void*) {
		vtkSmartPointer<vtkButtonWidget> buttonWidget = reinterpret_cast<vtkButtonWidget*>(caller);
		int state = buttonWidget->GetSliderRepresentation()->GetState();
		// actor on
		if (state == 0) {
			renderer->RemoveActor(actor);
			sliderWidget->EnabledOff();
		}
		// actor off
		else {
			renderer->AddActor(actor);
			sliderWidget->EnabledOn();
		}

	}
	vtkButtonCallback() : renderer(0), actor(0), renderWindow(0), sliderWidget(0) {}
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkSliderWidget> sliderWidget;
};


void readVTKFile()
{
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	/*
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\clams\\CLaMS_Nabro.vtk");
	reader->Update();
	*/
	

	
	vtkSmartPointer<vtkStructuredGridReader> earthReader = vtkSmartPointer<vtkStructuredGridReader>::New();
	earthReader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\support\\earth_texture.vtk");
	earthReader->Update();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> geometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	geometryFilter->SetInputConnection(earthReader->GetOutputPort());
	geometryFilter->Update();
	

	//vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
	//for (int i = 0; i < 10; i++)
	//	cout << data->GetPointData()->GetArrayName(i) << endl;
	//auto index = data->GetPointData()->GetArray("altitude");
	//data->GetPointData()->SetScalars(index);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(geometryFilter->GetOutput());
	//mapper->ScalarVisibilityOff();

	//vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	//mapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// clams
	/*
	vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
	reader2->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\clams\\CLaMS_Nabro.vtk");
	reader2->Update();

	vtkSmartPointer<vtkPolyData> data = reader2->GetOutput();

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetNumberOfTableValues(256);
	lut->Build();
	lut->SetAlpha(0.0);

	auto index = data->GetPointData()->GetArray("altitude");
	data->GetPointData()->SetScalars(index);

	vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputData(data);

	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->GetProperty()->SetOpacity(0.001);
	actor2->SetMapper(mapper2);
	*/




	// air_resampled
	vtkSmartPointer<vtkXMLImageDataReader> airReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	airReader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\airs_resampled\\uniform_ash.vti");
	airReader->Update();

	auto data = airReader->GetOutput();
	int *dim = data->GetDimensions();
	cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;

	vector<vtkSmartPointer<vtkPolyData>> airdata;

	vtkSmartPointer<vtkLookupTable> airlut = vtkSmartPointer<vtkLookupTable>::New();
	airlut->SetNumberOfTableValues(256);
	airlut->SetHueRange(0.0, 0.667);
	airlut->SetAlphaRange(0.0, 1.0);
	airlut->SetTableRange(-5, 0);
	airlut->Build();

	for (int i = 0; i < dim[2]; i++)
	{
		vtkSmartPointer<vtkExtractVOI> extractVOI = vtkSmartPointer<vtkExtractVOI>::New();
		extractVOI->SetInputData(data);
		extractVOI->SetVOI(0, dim[0] - 1, 0, dim[1] - 1, i, i);
		extractVOI->Update();

		auto extracted = extractVOI->GetOutput();
		vtkSmartPointer<vtkImageDataGeometryFilter> imageDataGeometryFilter = vtkSmartPointer<vtkImageDataGeometryFilter>::New();
		imageDataGeometryFilter->SetInputData(extracted);
		imageDataGeometryFilter->Update();
		vtkSmartPointer<vtkPolyData> data_tmp = imageDataGeometryFilter->GetOutput();
		auto index = data_tmp->GetPointData()->GetArray("ash");
		double range[2];
		index->GetRange(range, 0);
		cout << range[0] << " " << range[1] << endl;
		data_tmp->GetPointData()->SetScalars(index);
		for (int j = 0; j < data_tmp->GetNumberOfPoints(); j++)
		{
			double co[3];
			data_tmp->GetPoints()->GetPoint(j, co);
			data_tmp->GetPoints()->SetPoint(j, co[0], co[1], 0);
		}
		airdata.push_back(data_tmp);
	}

	vtkSmartPointer<vtkPolyDataMapper> airmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	airmapper->SetInputData(airdata[0]);
	airmapper->SetLookupTable(airlut);
	airmapper->UseLookupTableScalarRangeOn();

	vtkSmartPointer<vtkActor> airactor = vtkSmartPointer<vtkActor>::New();
	airactor->GetProperty()->SetOpacity(0.5);
	airactor->SetMapper(airmapper);





	// tropopause
	/*int troponum = 85;
	string troponame = "D:\\documents\\visualization\\project\\Volcanoes\\tropopause\\Tropoause-VISContest_000.vtk";
	vector<vtkSmartPointer<vtkPolyData>> tropodata;
	
	vtkSmartPointer<vtkLookupTable> tropolut = vtkSmartPointer<vtkLookupTable>::New();
	tropolut->SetNumberOfTableValues(256);
	tropolut->SetHueRange(0.0, 0.667);
	tropolut->SetAlphaRange(0.0, 1.0);
	tropolut->SetTableRange(5, 20);
	tropolut->Build();

	for (int i = 0; i < troponum; i++)
	{
		troponame[78] = '0' + i / 10;
		troponame[79] = '0' + i % 10;
		vtkSmartPointer<vtkStructuredGridReader> reader_tmp = vtkSmartPointer<vtkStructuredGridReader>::New();
		reader_tmp->SetFileName(troponame.c_str());
		reader_tmp->Update();
		vtkSmartPointer<vtkStructuredGridGeometryFilter> geometryFilter_tmp = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
		geometryFilter_tmp->SetInputConnection(reader_tmp->GetOutputPort());
		geometryFilter_tmp->Update();
		vtkSmartPointer<vtkPolyData> data = geometryFilter_tmp->GetOutput();
		auto index = data->GetPointData()->GetArray("trop_1");
		data->GetPointData()->SetScalars(index);
		tropodata.push_back(data);
	}

	vtkSmartPointer<vtkPolyDataMapper> tropomapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	tropomapper->SetInputData(tropodata[0]);
	tropomapper->SetLookupTable(tropolut);
	tropomapper->UseLookupTableScalarRangeOn();

	vtkSmartPointer<vtkActor> tropoactor = vtkSmartPointer<vtkActor>::New();
	tropoactor->GetProperty()->SetOpacity(1);
	tropoactor->SetMapper(tropomapper);*/



	// create a 2D slider
	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
	sliderRep->SetMinimumValue(0);
	//sliderRep->SetMaximumValue(84);
	sliderRep->SetMaximumValue(209);
	sliderRep->SetValue(0);
	//sliderRep->SetTitleText("tropo time");
	sliderRep->SetTitleText("air time");
	// set color properties
	sliderRep->GetSliderProperty()->SetColor(0.2, 0.2, 0.6);	// Change the color of the knob that slides
	sliderRep->GetTitleProperty()->SetColor(0, 0, 0);			// Change the color of the text indicating what the slider controls
	sliderRep->GetLabelProperty()->SetColor(0, 0, 0.4);			// Change the color of the text displaying the value
	sliderRep->GetSelectedProperty()->SetColor(0.4, 0.8, 0.4);	// Change the color of the knob when the mouse is held on it
	sliderRep->GetTubeProperty()->SetColor(0.7, 0.7, 0.7);		// Change the color of the bar
	sliderRep->GetCapProperty()->SetColor(0.7, 0.7, 0.7);		// Change the color of the ends of the bar
	// set position of the slider
	sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
	sliderRep->GetPoint1Coordinate()->SetValue(1000, 80);
	sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
	sliderRep->GetPoint2Coordinate()->SetValue(1160, 80);
	vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
	sliderWidget->SetInteractor(renWinInteractor);
	sliderWidget->SetRepresentation(sliderRep);
	sliderWidget->SetAnimationModeToAnimate();

	// create the callback
	vtkSmartPointer<vtkSliderCallback> callback = vtkSmartPointer<vtkSliderCallback>::New();
	//callback->mapper = tropomapper;
	//callback->data = tropodata;
	callback->mapper = airmapper;
	callback->data = airdata;
	// TODO: assign the marching cubes object (isosurface) to the marching cubes interaction callback
	sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);
	




	unsigned char green[3] = { 145,207,96 };
	unsigned char gray[3] = { 153,153,153 };
	vtkSmartPointer<vtkImageData> image1 = vtkSmartPointer<vtkImageData>::New();
	vtkSmartPointer<vtkImageData> image2 = vtkSmartPointer<vtkImageData>::New();
	CreateImage(image1, green, gray);
	CreateImage(image2, gray, green);

	// Create the widget and its representation
	vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
	buttonRepresentation->SetNumberOfStates(2);
	buttonRepresentation->SetButtonTexture(0, image1);
	buttonRepresentation->SetButtonTexture(1, image2);

	vtkSmartPointer<vtkButtonWidget> buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
	buttonWidget->SetInteractor(renWinInteractor);
	buttonWidget->SetRepresentation(buttonRepresentation);


	// Place the widget. Must be done after a render so that the viewport is defined.
	// Here the widget placement is in normalized display coordinates
	vtkSmartPointer<vtkCoordinate> upperLeft = vtkSmartPointer<vtkCoordinate>::New();
	upperLeft->SetCoordinateSystemToNormalizedDisplay();
	upperLeft->SetValue(0, 1.0);

	double bds[6];
	double sz = 50.0;
	bds[0] = upperLeft->GetComputedDisplayValue(renderer)[0] - sz;
	bds[1] = bds[0] + sz;
	bds[2] = upperLeft->GetComputedDisplayValue(renderer)[1] - sz;
	bds[3] = bds[2] + sz;
	bds[4] = bds[5] = 0.0;

	// Scale to 1, default is .5
	buttonRepresentation->SetPlaceFactor(1);
	buttonRepresentation->PlaceWidget(bds);
	buttonWidget->On();

	vtkSmartPointer<vtkButtonCallback> callbackButton = vtkSmartPointer<vtkButtonCallback>::New();
	callbackButton->renderer = renderer;
	callbackButton->actor = airactor;
	callbackButton->renderWindow = renWin;
	callbackButton->sliderWidget = sliderWidget;

	buttonWidget->AddObserver(vtkCommand::StateChangedEvent, callbackButton);
	sliderWidget->EnabledOn();

	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(0, 0, 1);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();

	renderer->AddActor(airactor);
	renderer->AddActor(actor);
	renderer->SetActiveCamera(aCamera);
	renderer->ResetCamera();
	renderer->SetBackground(.5, .5, .5);
	renderer->ResetCameraClippingRange();

	renWin->SetSize(1200, 800);
	renWin->Render();
	renWinInteractor->Initialize();
	renWinInteractor->Start();
}

int main()
{
	readVTKFile();
	return 0;
}