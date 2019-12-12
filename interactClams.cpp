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
#include <vtkCellData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkProperty.h>
#include <vtkExtractVOI.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkThresholdPoints.h>
#include <vtkExtractSelection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkInformation.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkRegularPolygonSource.h>
#include <vtkGlyph2D.h>

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <set>

using namespace std;

const int time_interval = 86400;

void getColorCorrespondingTovalue(double val, double& r, double& g, double& b, double range)
{
	static const int numColorNodes = 9;
	double color[numColorNodes][3] =
	{
		255 / 255.0,255 / 255.0,255 / 255.0,
		247 / 255.0,235 / 255.0,222 / 255.0,
		239 / 255.0,219 / 255.0,198 / 255.0,
		225 / 255.0,202 / 255.0,158 / 255.0,
		214 / 255.0,174 / 255.0,107 / 255.0,
		198 / 255.0,146 / 255.0,66 / 255.0,
		181 / 255.0,113 / 255.0,33 / 255.0,
		156 / 255.0,81 / 255.0,8 / 255.0,
		107 / 255.0,48 / 255.0,8 / 255.0
	};

	for (int i = 0; i < (numColorNodes - 1); i++)
	{
		double currFloor = 0 + ((double)i / (numColorNodes - 1)) * range;
		double currCeil = 0 + ((double)(i + 1) / (numColorNodes - 1)) * range;

		if ((val >= currFloor) && (val <= currCeil))
		{
			double currFraction = (val - currFloor) / (currCeil - currFloor);
			r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
			g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
			b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
		}
	}
}

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


class vtkSliderCallbackclams : public vtkCommand
{
public:
	static vtkSliderCallbackclams* New() {
		return new vtkSliderCallbackclams;
	}
	virtual void Execute(vtkObject* caller, unsigned long, void*) {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		double now = static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue();

		auto times = data->GetPointData()->GetArray("time");
		double timerange[2];
		times->GetRange(timerange, 0);
		now += timerange[0];

		auto altitude = data->GetPointData()->GetArray("altitude");

		vtkSmartPointer<vtkFloatArray> newindex = vtkSmartPointer<vtkFloatArray>::New();
		newindex->SetNumberOfComponents(1);
		for (int i = 0; i < data->GetNumberOfPoints(); i++)
		{
			double time = times->GetTuple1(i);
			if (time >= now && time < now + time_interval)
			{
				newindex->InsertNextValue(altitude->GetTuple1(i) + 2);
			}
			else
			{
				newindex->InsertNextValue(0);
			}
		}
		data->GetPointData()->SetScalars(newindex);
	}
	vtkSliderCallbackclams() : data(0) {}
	vtkSmartPointer<vtkPolyData> data;
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
		// actor plane
		if (state == 0) {
			mapper->SetInputData(planedata);
			emapper->SetInputData(eplanedata);
			slider->data = planedata;
		}
		// actor sphere
		else {
			mapper->SetInputData(spheredata);
			emapper->SetInputData(espheredata);
			slider->data = spheredata;
		}

	}
	vtkButtonCallback() : planedata(0), spheredata(0), mapper(0), eplanedata(0), espheredata(0), emapper(0), slider(0) {}
	vtkSmartPointer<vtkPolyData> planedata;
	vtkSmartPointer<vtkPolyData> spheredata;
	vtkSmartPointer<vtkPolyData> eplanedata;
	vtkSmartPointer<vtkPolyData> espheredata;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkPolyDataMapper> emapper;
	vtkSmartPointer<vtkSliderCallbackclams> slider;
};


void carte2sphere(vtkSmartPointer<vtkPolyData> data, double r)
{
	double pi = 3.1415926535897;
	for (int i = 0; i < data->GetNumberOfPoints(); i++)
	{
		double co[3];
		data->GetPoints()->GetPoint(i, co);
		double z = r * sin(co[1] / 180 * pi);
		double x = r * cos(co[1] / 180 * pi) * cos(co[0] / 180 * pi);
		double y = r * cos(co[1] / 180 * pi) * sin(co[0] / 180 * pi);
		data->GetPoints()->SetPoint(i, x, y, z);
	}
}


void readVTKFile()
{
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->PointSmoothingOn();
	cout << renWin->GetPointSmoothing() << endl;
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(0, 0, 1);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();



	// background
	vtkSmartPointer<vtkStructuredGridReader> earthReader = vtkSmartPointer<vtkStructuredGridReader>::New();
	earthReader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\support\\earth_texture.vtk");
	earthReader->Update();
	vtkSmartPointer<vtkStructuredGridGeometryFilter> geometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	geometryFilter->SetInputConnection(earthReader->GetOutputPort());
	geometryFilter->Update();

	vtkSmartPointer<vtkPolyData> earthdata = geometryFilter->GetOutput();
	for (int i = 0; i < earthdata->GetNumberOfPoints(); i++)
	{
		double co[3];
		earthdata->GetPoints()->GetPoint(i, co);
		earthdata->GetPoints()->SetPoint(i, co[0], co[1], -0.1);
	}
	vtkSmartPointer<vtkPolyData> earthdata_ori = vtkSmartPointer<vtkPolyData>::New();
	earthdata_ori->DeepCopy(earthdata);
	carte2sphere(earthdata, 100);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(earthdata);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// clams

	vtkSmartPointer<vtkPolyDataReader> clamsreader = vtkSmartPointer<vtkPolyDataReader>::New();
	clamsreader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\clams\\CLaMS_Puyehue.vtk");
	//clamsreader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\clams\\CLaMS_Nabro.vtk");
	clamsreader->Update();

	vtkSmartPointer<vtkPolyData> clamsdata = clamsreader->GetOutput();

	auto trajs = clamsdata->GetCellData()->GetArray("seed_id");
	double trajrange[2];
	trajs->GetRange(trajrange, 0);
	cout << trajrange[0] << " " << trajrange[1] << endl;

	auto times = clamsdata->GetPointData()->GetArray("time");
	double timerange[2];
	times->GetRange(timerange, 0);
	cout << timerange[0] << " " << timerange[1] << endl;

	auto altitude = clamsdata->GetPointData()->GetArray("altitude");
	double altirange[2];
	altitude->GetRange(altirange, 0);
	cout << altirange[0] << " " << altirange[1] << endl;

	vtkSmartPointer<vtkFloatArray> newindex = vtkSmartPointer<vtkFloatArray>::New();
	newindex->SetNumberOfComponents(1);
	for (int i = 0; i < clamsdata->GetNumberOfPoints(); i++)
	{
		double time = times->GetTuple1(i);
		if (time < timerange[0] + time_interval)
		{
			newindex->InsertNextValue(altitude->GetTuple1(i) + 2);
		}
		else
		{
			newindex->InsertNextValue(0);
		}
	}
	clamsdata->GetPointData()->SetScalars(newindex);

	double range = altirange[1] + 2;
	int numColors = 256;
	vtkSmartPointer<vtkLookupTable> clamslut = vtkSmartPointer<vtkLookupTable>::New();
	clamslut->SetNumberOfTableValues(256);
	clamslut->SetScaleToLinear();
	double r, g, b;
	for (int i = 0; i < numColors; i++)
	{
		double val = 0 + ((double)i / numColors) * range;
		getColorCorrespondingTovalue(val, r, g, b, range);
		if (r == 1 && g == 1 && b == 1) { clamslut->SetTableValue(i, r, g, b, 0.0); }
		else { clamslut->SetTableValue(i, r, g, b, 1.0); }
	}
	cout << clamslut->GetOpacity(0) << " " << clamslut->GetOpacity(altirange[0] + 2) << endl;

	vtkSmartPointer<vtkPolyData> clamsdata_ori = vtkSmartPointer<vtkPolyData>::New();
	clamsdata_ori->DeepCopy(clamsdata);
	carte2sphere(clamsdata, 101);

	vtkSmartPointer<vtkPolyDataMapper> clamsmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	clamsmapper->SetInputData(clamsdata);
	clamsmapper->SetLookupTable(clamslut);
	clamsmapper->UseLookupTableScalarRangeOff();
	clamsmapper->SetScalarRange(0, altirange[1] + 2);

	vtkSmartPointer<vtkActor> clamsactor = vtkSmartPointer<vtkActor>::New();
	clamsactor->GetProperty()->SetOpacity(0.5);
	clamsactor->SetMapper(clamsmapper);



	// create a 2D slider
	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
	sliderRep->SetMinimumValue(0);
	sliderRep->SetMaximumValue(timerange[1] - timerange[0] - time_interval);
	sliderRep->SetValue(0);
	sliderRep->SetTitleText("time");
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
	vtkSmartPointer<vtkSliderCallbackclams> callback = vtkSmartPointer<vtkSliderCallbackclams>::New();
	callback->data = clamsdata;
	sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);




	// button
	unsigned char green[3] = { 145,207,96 };
	unsigned char gray[3] = { 153,153,153 };
	vtkSmartPointer<vtkImageData> image1 = vtkSmartPointer<vtkImageData>::New();
	vtkSmartPointer<vtkImageData> image2 = vtkSmartPointer<vtkImageData>::New();
	CreateImage(image1, green, gray);
	CreateImage(image2, gray, green);

	vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation = vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
	buttonRepresentation->SetNumberOfStates(2);
	buttonRepresentation->SetButtonTexture(0, image1);
	buttonRepresentation->SetButtonTexture(1, image2);

	vtkSmartPointer<vtkButtonWidget> buttonWidget = vtkSmartPointer<vtkButtonWidget>::New();
	buttonWidget->SetInteractor(renWinInteractor);
	buttonWidget->SetRepresentation(buttonRepresentation);

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



	// set callback
	vtkSmartPointer<vtkButtonCallback> callbackButton = vtkSmartPointer<vtkButtonCallback>::New();
	callbackButton->mapper = clamsmapper;
	callbackButton->planedata = clamsdata_ori;
	callbackButton->spheredata = clamsdata;
	callbackButton->emapper = mapper;
	callbackButton->eplanedata = earthdata_ori;
	callbackButton->espheredata = earthdata;
	callbackButton->slider = callback;

	buttonWidget->AddObserver(vtkCommand::StateChangedEvent, callbackButton);
	sliderWidget->EnabledOn();




	// general
	renderer->AddActor(clamsactor);
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
