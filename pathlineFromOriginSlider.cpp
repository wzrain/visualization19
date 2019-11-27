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

using namespace std;

const int time_interval = 86400;
static const long timeStart = 360460800; // P
static const long timeEnd = 367664384; // P
// static const long timeStart = 361217800; // N
// static const long timeEnd = 368582384; // N
static const double X = -72.0f; // P
static const double Y = -40.0f; // P
// static const double X = 41.0f; // N
// static const double Y = 13.0f; // N
static const double delta = 2.0f;
static const double timeDelta = 6 * 3600;
static const long slideGap = 6*3600;
static const long timeGap = 3600;
static const long timemax = 30*24*3600;
std::vector<bool> pointBool(1447514, false);

double timeStartHour(double h){
	return timeStart + h * 3600;
}

void getColorCorrespondingTovalue(double val, double& r, double& g, double& b, double range)
{
	static const int numColorNodes = 9;
	double color[numColorNodes][3] =
	{
		// 255 / 255.0,255 / 255.0,255 / 255.0,
		// 247 / 255.0,235 / 255.0,222 / 255.0,
		// 239 / 255.0,219 / 255.0,198 / 255.0,
		// 225 / 255.0,202 / 255.0,158 / 255.0,
		// 214 / 255.0,174 / 255.0,107 / 255.0,
		// 198 / 255.0,146 / 255.0,66 / 255.0,
		// 181 / 255.0,113 / 255.0,33 / 255.0,
		// 156 / 255.0,81 / 255.0,8 / 255.0,
		// 107 / 255.0,48 / 255.0,8 / 255.0

		255 / 255.0,255 / 255.0,255 / 255.0,
		254 / 255.0,224 / 255.0,210 / 255.0,
		252 / 255.0,187 / 255.0,161 / 255.0,
		252 / 255.0,146 / 255.0,114 / 255.0,
		251 / 255.0,106 / 255.0,74 / 255.0,
		239 / 255.0,59 / 255.0,44 / 255.0,
		203 / 255.0,24 / 255.0,29 / 255.0,
		165 / 255.0,15 / 255.0,21 / 255.0,
		103 / 255.0,0 / 255.0,13 / 255.0,
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


class vtkSliderCallbackclams : public vtkCommand
{
public:
	static vtkSliderCallbackclams* New() {
		return new vtkSliderCallbackclams;
	}
	virtual void Execute(vtkObject* caller, unsigned long, void*) {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		double now = static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue();
		long longnow = long(now) / slideGap * timeGap;

		auto times = data->GetPointData()->GetArray("time");
		double timerange[2];
		times->GetRange(timerange, 0);
		longnow += timerange[0];

		auto altitude = data->GetPointData()->GetArray("altitude");

		vtkSmartPointer<vtkFloatArray> newindex = vtkSmartPointer<vtkFloatArray>::New();
		newindex->SetNumberOfComponents(1);
		for (int i = 0; i < data->GetNumberOfPoints(); i++){
			double time = times->GetTuple1(i);
			if (time < longnow and pointBool[i]){
				newindex->InsertNextValue(altitude->GetTuple1(i) + 2);
			}
			else{
				newindex->InsertNextValue(0);
			}
		}
		data->GetPointData()->SetScalars(newindex);
	}
	vtkSliderCallbackclams() : data(0) {}
	vtkSmartPointer<vtkPolyData> data;
};


class vtkSliderCallbackmipas : public vtkCommand
{
public:
	static vtkSliderCallbackmipas* New() {
		return new vtkSliderCallbackmipas;
	}
	virtual void Execute(vtkObject* caller, unsigned long, void*) {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		double now = static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue();

		auto mipastimes = data->GetPointData()->GetArray("time");
		double mipastrange[2];
		mipastimes->GetRange(mipastrange, 0);
		now += mipastrange[0];

		vtkSmartPointer<vtkThresholdPoints> mipasthres = vtkSmartPointer<vtkThresholdPoints>::New();
		mipasthres->SetInputData(data);
		mipasthres->ThresholdBetween(1.9, 2.1);
		mipasthres->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "detection");
		mipasthres->Update();
		auto mipasash = mipasthres->GetOutput();

		vtkSmartPointer<vtkThresholdPoints> mipasthres2 = vtkSmartPointer<vtkThresholdPoints>::New();
		mipasthres2->SetInputData(mipasash);
		mipasthres2->ThresholdBetween(now, now + time_interval);
		mipasthres2->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "time");
		mipasthres2->Update();
		auto filtered = mipasthres2->GetOutput();

		vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
		polygonSource->SetNumberOfSides(10);
		polygonSource->SetRadius(1);
		vtkSmartPointer<vtkGlyph2D> glyph2D = vtkSmartPointer<vtkGlyph2D>::New();
		glyph2D->SetSourceConnection(polygonSource->GetOutputPort());
		glyph2D->SetInputData(filtered);
		glyph2D->Update();

		mapper->SetInputData(glyph2D->GetOutput());
	}
	vtkSliderCallbackmipas() : data(0) {}
	vtkSmartPointer<vtkPolyData> data;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
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
	renWin->PointSmoothingOn();
	renWin->GetPointSmoothing();
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	// background
	vtkSmartPointer<vtkStructuredGridReader> earthReader = vtkSmartPointer<vtkStructuredGridReader>::New();
	earthReader->SetFileName("../../../support/earth_texture.vtk");
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

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(earthdata);
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// clams	
	vtkSmartPointer<vtkPolyDataReader> clamsreader = vtkSmartPointer<vtkPolyDataReader>::New();
	clamsreader->SetFileName("../../../clams/CLaMS_Puyehue.vtk");
	clamsreader->Update();
	vtkSmartPointer<vtkPolyData> clamsdata = clamsreader->GetOutput();

	auto trajs = clamsdata->GetCellData()->GetArray("seed_id");
	double trajrange[2];
	trajs->GetRange(trajrange, 0);

	auto times = clamsdata->GetPointData()->GetArray("time");
	double timerange[2];
	times->GetRange(timerange, 0);

	auto altitude = clamsdata->GetPointData()->GetArray("altitude");
	double altirange[2];
	altitude->GetRange(altirange, 0);

	vtkSmartPointer<vtkFloatArray> newindex = vtkSmartPointer<vtkFloatArray>::New();
	newindex->SetNumberOfValues(clamsdata->GetNumberOfPoints());
	for(int i = 0; i < clamsdata->GetNumberOfPoints(); i++){newindex->InsertNextValue(0);}
	for(int i = 0; i < trajs->GetNumberOfTuples(); i++){
		int pointNum = clamsdata->GetCell(i)->GetPointIds()->GetNumberOfIds();
		std::vector<vtkIdType> PIdx;
		std::vector<double> PX;
		bool near = false;
		for(int j = 0; j < pointNum; j++){
			double pTime = times->GetTuple1(clamsdata->GetCell(i)->GetPointId(j)); 
			double co[3];
			clamsdata->GetPoints()->GetPoint(clamsdata->GetCell(i)->GetPointId(j),co);
			if(!near and std::abs(co[0]-X) < delta and std::abs(co[1]-Y) < delta and pTime < timerange[0] + timeDelta){near = true;}
			PIdx.push_back(clamsdata->GetCell(i)->GetPointId(j));
			PX.push_back(co[0]);
		}
		if(near){
			for(int j = 0; j < PX.size(); j++){
				if(PX[j] > -170 and PX[j] < 170){
					pointBool[PIdx[j]] = true;
				}
			}
		}
	}
	clamsdata->GetPointData()->SetScalars(newindex);

	// look up table
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

	// mapper
	vtkSmartPointer<vtkPolyDataMapper> clamsmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	clamsmapper->SetInputData(clamsdata);
	clamsmapper->SetLookupTable(clamslut);
	clamsmapper->UseLookupTableScalarRangeOff();
	clamsmapper->SetScalarRange(0, altirange[1] + 2);

	// actor
	vtkSmartPointer<vtkActor> clamsactor = vtkSmartPointer<vtkActor>::New();
	clamsactor->GetProperty()->SetOpacity(0.5);
	clamsactor->SetMapper(clamsmapper);
		
	// create a 2D slider
	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
	sliderRep->SetMinimumValue(0);
	//sliderRep->SetMaximumValue(84);
	//sliderRep->SetMaximumValue(209);
	sliderRep->SetMaximumValue(timerange[1] - timerange[0] - time_interval);
	// sliderRep->SetMaximumValue(mipastrange[1] - mipastrange[0] - time_interval);
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
	sliderWidget->EnabledOn();
	
	// general
	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(0, 0, 1);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();

	//renderer->AddActor(tropoactor);
	//renderer->AddActor(airactor);
	renderer->AddActor(clamsactor);
	// renderer->AddActor(mipasactor);
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