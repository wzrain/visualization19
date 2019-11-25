#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCoordinate.h>
#include <vtkSphereSource.h>
#include <vtkButtonWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridGeometryFilter.h>
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
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <iostream>
#include <vector>
#include <string>

//https://math.stackexchange.com/questions/915290/find-vector-field-given-curl

static const double min = 0;
static const double max = 650;
static const double range = max-min;
static const double numColors = 100;
static const long timeStart = 360460800;
static const long timeEnd = 367664384;

void getColorCorrespondingTovalue(double val, double &r, double &g, double &b)
{
	static const int numColorNodes = 9;
	double color[numColorNodes][3] =
	{
		255/255.0,255/255.0,255/255.0,
		222/255.0,235/255.0,247/255.0,
		198/255.0,219/255.0,239/255.0,
		158/255.0,202/255.0,225/255.0,
		107/255.0,174/255.0,214/255.0,
		66/255.0,146/255.0,198/255.0,
		33/255.0,113/255.0,181/255.0,
		8/255.0,81/255.0,156/255.0,
		8/255.0,48/255.0,107/255.0
	};

	for (int i = 0; i < (numColorNodes - 1); i++)
	{
		double currFloor = min + ((double)i / (numColorNodes - 1)) * range;
		double currCeil = min + ((double)(i + 1) / (numColorNodes - 1)) * range;

		if ((val >= currFloor) && (val <= currCeil))
		{
			double currFraction = (val - currFloor) / (currCeil - currFloor);
			r = color[i][0] * (1.0 - currFraction) + color[i + 1][0] * currFraction;
			g = color[i][1] * (1.0 - currFraction) + color[i + 1][1] * currFraction;
			b = color[i][2] * (1.0 - currFraction) + color[i + 1][2] * currFraction;
		}
	}
}

long timeStartHour(long h){
	return timeStart + h * 3600;
}

void readVTKFile(const char *fileName,const char *attribute)
{
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	if(std::string(fileName) == std::string("")){
		// reader->SetFileName("../../../clams/CLaMS_Nabro.vtk");
		reader->SetFileName("../../../clams/CLaMS_Puyehue.vtk");
	}
	else{reader->SetFileName(fileName);}
	reader->Update();

	vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
	// vtkSmartPointer<vtkUnstructuredGrid> ugrid = reader->GetOutput();
    // vtkSmartPointer<vtkCellData> cellData = data->GetCellData();
    // vtkSmartPointer<vtkDataArray> celldata = cellData->GetScalars();
	// vtkSmartPointer<vtkCellArray> cell = reader->GetOutput();
	// cout << data->GetPointData()->GetArrayName(0) << endl;
	// cout << data->GetPointData()->GetArrayName(1) << endl;
	// cout << data->GetPointData()->GetArrayName(2) << endl;
	// cout << data->GetPointData()->GetArrayName(3) << endl;
	// cout << data->GetPointData()->GetArrayName(4) << endl;
	// cout << data->GetPointData()->GetArrayName(5) << endl;
	std::cout << data->GetNumberOfCells() << std::endl;
	std::cout << data->GetNumberOfPoints() << std::endl;
	vtkSmartPointer<vtkFloatArray> idx = vtkSmartPointer<vtkFloatArray>::New();
	idx->SetNumberOfValues(data->GetNumberOfPoints());
	for(int i = 0; i < data->GetNumberOfPoints(); i++){
		idx->SetValue(i,0);
	}
	// std::cout << data->GetCellData()->GetNumberOfCells() << std::endl;
	// auto trajs = data->GetCellData()->GetArray("seed_id");
	// data->GetCellData()->SetScalars(trajs);
	// for(int i = 0; i < data->GetNumberOfCells(); i++){
	// 	std::cout << i << ' ' << celldata->GetTuple1(i) << std::endl;
	// }
	if(std::string(attribute) == std::string("")){
		auto time = data->GetPointData()->GetArray("time");
		auto altitude = data->GetPointData()->GetArray("altitude");
		auto pressure = data->GetPointData()->GetArray("pressure");
		auto temperature = data->GetPointData()->GetArray("temperature");
		auto potTemperature = data->GetPointData()->GetArray("pot_temperature");
		auto potVorticity = data->GetPointData()->GetArray("pot_vorticity");
		auto trajs = data->GetCellData()->GetArray("seed_id"); // the cell type is VTK_POLY_LINE
		vtkSmartPointer<vtkCellArray> trajsThreshold = vtkSmartPointer<vtkCellArray>::New();
		// std::vector<std::vector<int>> times(trajs->GetNumberOfTuples());
		std::cout << trajs->GetNumberOfTuples() << std::endl;
		// long timeStart = long(1 << 28) * long(100);
		// long timeEnd = 0;
		double co[3];
		data->GetPoints()->GetPoint(1,co);
		long timeThresholdUp = timeStartHour(24*1);
		long timeThresholdDown = timeStartHour(0);
		for(int i = 0; i < trajs->GetNumberOfTuples(); i++){
			int pointNum = data->GetCell(i)->GetPointIds()->GetNumberOfIds();
			// std::cout << pointNum << std::endl;
			bool add = false;
			int addCnt = 0;
			std::vector<double> PTime;
			std::vector<vtkIdType> PIdx;
			std::vector<double> PY;
			for(int j = 0; j < pointNum; j++){
				float pTime = time->GetTuple1(data->GetCell(i)->GetPointId(j));
				if(pTime < timeThresholdUp and pTime > timeThresholdDown){
					// double co[3];
					// data->GetPoints()->GetPoint(data->GetCell(i)->GetPointId(j),co);
					// PTime.push_back(trajs->GetTuple1(i));
					// PIdx.push_back(data->GetCell(i)->GetPointId(j));
					// if(!add and PY.size()){
					// 	if(std::abs(co[1] - PY.back()) > 3){add = true;}
					// }
					// PY.push_back(co[1]);
					idx->SetValue(data->GetCell(i)->GetPointId(j),trajs->GetTuple1(i));
				}
				// if(pTime == 0){std::cout << i << ' ' << j << std::endl;}
				// times[i].push_back(pTime);
				// if(timeStart > pTime){timeStart = pTime;}
				// if(timeEnd < pTime){timeEnd = pTime;}
			}
			// if(add){
			// 	for(int j = 0; j < PTime.size(); j++){
			// 		idx->SetValue(PIdx[j],PTime[j]);
			// 	}
			// }
			// std::cout << i << ' ' << trajs->GetTuple1(i) << std::endl;
			// std::cout << i << ' ' << data->GetCell(i)->GetCellType() << std::endl;
			// // std::cout << i << ' ' << data->GetCell(i)->GetPoints()->GetNumberOfPoints() << std::endl;
			// std::cout << i << ' ' << data->GetCell(i)->GetPointIds()->GetNumberOfIds() << std::endl;
			// std::cout << i << ' ' << data->GetCell(i)->GetPointId(0) << std::endl;
			// std::cout << i << ' ' << time->GetTuple1(data->GetCell(i)->GetPointId(0))<< std::endl;
			// std::cout << "================\n";
		}
		// std::cout << long(timeStart) << ' ' << long(timeEnd) << std::endl;
		// data->GetPointData()->SetScalars(time);
		// data->GetPointData()->SetScalars(altitude);
		// data->GetPointData()->SetScalars(pressure);
		// data->GetPointData()->SetScalars(temperature);
		// data->GetPointData()->SetScalars(potTemperature);
		// data->GetPointData()->SetScalars(potVorticity);
		data->GetPointData()->SetScalars(idx);
		// data->GetCellData()->SetScalars(trajs);
	}
	else{
		auto ashIndex = data->GetPointData()->GetArray(attribute);
		data->GetPointData()->SetScalars(ashIndex);
	}
	
	vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	lookupTable->SetScaleToLinear();
	lookupTable->SetNumberOfTableValues(numColors);
	double r, g, b;
	for (int i = 0; i < numColors; i++){
		double val = min + ((double)i / numColors) * range;
		getColorCorrespondingTovalue(val, r, g, b);
		if(r == 1 and g == 1 and b == 1){lookupTable->SetTableValue(i, r, g, b, 0.0);}
		else{lookupTable->SetTableValue(i, r, g, b, 1.0);}
	}	

	vtkSmartPointer <vtkScalarBarActor> legend = vtkScalarBarActor::New();
	legend->SetLookupTable(lookupTable);
	legend->SetNumberOfLabels(5);
	legend->SetTitle("id of trajectory");
	legend->SetVerticalTitleSeparation(6);
	legend->GetPositionCoordinate()->SetValue(0.88,0.1);
	legend->SetWidth(0.05);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(reader->GetOutputPort());
	mapper->SetLookupTable(lookupTable);
	mapper->SetScalarRange(0,650);
	//mapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
	aCamera->SetViewUp(0, 1, 0);
	aCamera->SetPosition(0, 0, 1);
	aCamera->SetFocalPoint(0, 0, 0);
	aCamera->ComputeViewPlaneNormal();

	renderer->AddActor(actor);
	renderer->AddActor2D(legend);
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
