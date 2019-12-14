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
#include <vtkVertexGlyphFilter.h>
#include <vtkScalarBarActor.h>

#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <map>
// 361224000
// 361245600

using namespace std;

struct Node {
	vtkIdType pid, tid;
	double alt;
	bool above;
	Node() {}
	Node(vtkIdType pi, vtkIdType ti, double t, bool flag) { pid = pi; tid = ti; alt = t; above = flag; }
};

const int time_interval = 86400;
// static const long timeStart = 360460800; // P
// static const long timeEnd = 367664384; // P
// static const long timeStart = 361217800; // N
// static const long timeEnd = 368582384; // N
// static const double X = -72.0f; // P
// static const double Y = -40.0f; // P
static const double NX = 41.0f; // N
static const double NY = 13.0f; // N
static const double GX = -17.0f; // G
static const double GY = 64.0f; // G
static const double delta = 3.0f;
static const double timeDelta = 6 * 3600;
static const long slideGap = 5 * 3600;
static const long timeGap = 3600;
static const long timemax = 30 * 24 * 3600;
static const long timePath = 6 * 3600;
double visualMin;
double visualMax;
long timeStart, timeEnd;
#define visualRange (visualMax-visualMin)
// std::vector<bool> pointBool(1447514, false);
std::vector<int> pointBool(24000661, 0); // G
std::map<long, std::vector<Node>> timeLine;
std::vector<long> trajStartTime;

double timeStartHour(double h) {
	return timeStart + h * 3600;
}

void getColorCorrespondingTovalue(double val, double& r, double& g, double& b, double range)
{
	static const int numColorNodes = 9;
	double color[numColorNodes][3] =
	{
		0 / 255.0,0 / 255.0,255 / 255.0,
		64 / 255.0,64 / 255.0,192 / 255.0,
		128 / 255.0,128 / 255.0,128 / 255.0,
		192 / 255.0,192 / 255.0,64 / 255.0,
		255 / 255.0,255 / 255.0,0 / 255.0,
		255 / 255.0,192 / 255.0,0 / 255.0,
		255 / 255.0,128 / 255.0,0 / 255.0,
		255 / 255.0,64 / 255.0,0 / 255.0,
		255 / 255.0,0 / 255.0,0 / 255.0,
	};

	for (int i = 0; i < (numColorNodes - 1); i++)
	{
		double currFloor = visualMin + ((double)i / (numColorNodes - 1)) * visualRange;
		double currCeil = visualMin + ((double)(i + 1) / (numColorNodes - 1)) * visualRange;

		if ((val >= currFloor) && (val <= currCeil))
		{
			double currFraction = (val - currFloor) / (currCeil - currFloor);
			r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
			g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
			b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
		}
	}
}


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

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(earthdata);
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	// clams	
	vtkSmartPointer<vtkPolyDataReader> clamsreader = vtkSmartPointer<vtkPolyDataReader>::New();
	clamsreader->SetFileName("D:\\documents\\visualization\\project\\Volcanoes\\clams\\CLaMS_Nabro.vtk");
	clamsreader->Update();
	vtkSmartPointer<vtkPolyData> clamsdata = clamsreader->GetOutput();

	auto trajs = clamsdata->GetCellData()->GetArray("seed_id");
	double trajrange[2];
	trajs->GetRange(trajrange, 0);

	auto times = clamsdata->GetPointData()->GetArray("time");
	double timerange[2];
	times->GetRange(timerange, 0);
	timeStart = timerange[0]; timeEnd = timerange[1];


	// consider other feature
	map<string, pair<double, double>> featminmax;
	featminmax["pressure"] = make_pair(1.6, 7.2);
	featminmax["temperature"] = make_pair(1.1, 2.2);
	featminmax["pot_vorticity"] = make_pair(0.5, 0.9);
	featminmax["pot_temperature"] = make_pair(0.1, 0.5);

	string featname = "temperature";
	auto clamsfeat = clamsdata->GetPointData()->GetArray(featname.c_str());
	double featrange[2];
	clamsfeat->GetRange(featrange, 0);
	cout << featrange[0] << " " << featrange[1] << endl;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> gradindex = vtkSmartPointer<vtkFloatArray>::New();
	double globalmax = 0;
	for (int i = 0; i < clamsdata->GetNumberOfCells(); i++)
	{
		if (i % 5000 == 0) { std::cout << i << std::endl; }
		int numpoints = clamsdata->GetCell(i)->GetPointIds()->GetNumberOfIds();
		double maxgrad = 0;
		int maxpid = clamsdata->GetCell(i)->GetPointId(0);
		bool near = false;
		for (int j = 1; j < numpoints - 1; j++)
		{
			int pointid = clamsdata->GetCell(i)->GetPointId(j);
			int pointid1 = clamsdata->GetCell(i)->GetPointId(j - 1);
			int pointid2 = clamsdata->GetCell(i)->GetPointId(j + 1);
			double featdelta = clamsfeat->GetTuple1(pointid2) - clamsfeat->GetTuple1(pointid1);
			double timedelta = times->GetTuple1(pointid2) - times->GetTuple1(pointid1);
			double grad = fabs(featdelta) * 3600 / fabs(timedelta);
			if (grad > maxgrad)
			{
				maxgrad = grad;
				maxpid = pointid;
			}
			double cotmp[3];
			clamsdata->GetPoint(pointid, cotmp);
			if (std::abs(cotmp[0] - NX) <= delta && std::abs(cotmp[1] - NY) <= delta)
			{
				near = true;
			}
		}

		if (near && maxgrad > featminmax[featname].first)
		{
			double co[3];
			clamsdata->GetPoint(maxpid, co);
			points->InsertNextPoint(co);
			gradindex->InsertNextTuple1(maxgrad);
			if (maxgrad > globalmax)
			{
				globalmax = maxgrad;
			}
		}
	}
	cout << globalmax << endl;

	vtkSmartPointer<vtkPolyData> pointpoly = vtkSmartPointer<vtkPolyData>::New();
	pointpoly->SetPoints(points);
	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointpoly);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> newdata = vtkSmartPointer<vtkPolyData>::New();
	newdata->ShallowCopy(vertexFilter->GetOutput());
	newdata->GetPointData()->SetScalars(gradindex);

	cout << newdata->GetNumberOfPoints() << endl;

	// look up table
	double range = globalmax + 2;
	int numColors = 256;
	vtkSmartPointer<vtkLookupTable> clamslut = vtkSmartPointer<vtkLookupTable>::New();
	clamslut->SetNumberOfTableValues(256);
	clamslut->SetScaleToLinear();
	double r, g, b;
	visualMin = featminmax[featname].first;
	visualMax = featminmax[featname].second;
	for (int i = 0; i < numColors; i++)
	{
		double val = visualMin + ((double)i / numColors) * visualRange;
		getColorCorrespondingTovalue(val, r, g, b, range);
		if (r == 1 && g == 1 && b == 1) { clamslut->SetTableValue(i, r, g, b, 0.0); }
		else { clamslut->SetTableValue(i, r, g, b, (r + 1) / 2); }
	}

	vtkSmartPointer<vtkPolyDataMapper> clamsnewmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	clamsnewmapper->SetInputData(newdata);
	clamsnewmapper->SetLookupTable(clamslut);
	clamsnewmapper->SetScalarRange(visualMin, visualMax);

	// actor
	vtkSmartPointer<vtkActor> clamsnewactor = vtkSmartPointer<vtkActor>::New();
	clamsnewactor->SetMapper(clamsnewmapper);
	clamsnewactor->GetProperty()->SetPointSize(2);

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(clamsnewmapper->GetLookupTable());
	scalarBar->SetNumberOfLabels(4);





	// general
	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(90, 45, 1);
	aCamera->SetFocalPoint(90, 45, 0);
	aCamera->ComputeViewPlaneNormal();

	//renderer->AddActor(tropoactor);
	//renderer->AddActor(airactor);
	//renderer->AddActor(clamsactor);
	renderer->AddActor(clamsnewactor);
	// renderer->AddActor(mipasactor);
	renderer->AddActor(actor);
	renderer->AddActor2D(scalarBar);
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