#include <string>
#include <vector>
#include <unordered_map>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkCoordinate.h>
#include <vtkSphereSource.h>
#include <vtkButtonWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkFloatArray.h>
#include <vtkProperty.h>
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
#include <vtkCellData.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkLookupTable.h>
#include <vtkDelaunay2D.h>
#include <vtkPolyLine.h>

static const double min = 0;
static const double max = 650;
static const double range = max - min;
static const double numColors = 100;
static const int nabro_x = 41;
static const int nabro_y = 13;

void getColorCorrespondingTovalue(double val, double &r, double &g, double &b)
{
	static const int numColorNodes = 9;
	double color[numColorNodes][3] =
	{
		//0.6980, 0.0941, 0.1686, // Red
		//0.8392, 0.3765, 0.3020,
		//0.9569, 0.6471, 0.5098,
		//0.9922, 0.8588, 0.7804,
		//0.9686, 0.9686, 0.9686, // White
		//0.8196, 0.8980, 0.9412,
		//0.5725, 0.7725, 0.8706,
		//0.2627, 0.5765, 0.7647,
		//0.1294, 0.4000, 0.6745  // Blue
		255 / 255.0,255 / 255.0,255 / 255.0,
		222 / 255.0,235 / 255.0,247 / 255.0,
		198 / 255.0,219 / 255.0,239 / 255.0,
		158 / 255.0,202 / 255.0,225 / 255.0,
		107 / 255.0,174 / 255.0,214 / 255.0,
		66 / 255.0,146 / 255.0,198 / 255.0,
		33 / 255.0,113 / 255.0,181 / 255.0,
		8 / 255.0,81 / 255.0,156 / 255.0,
		8 / 255.0,48 / 255.0,107 / 255.0
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

void readVTKFile()
{
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renWin);

	vtkSmartPointer<vtkPolyDataReader> clams_reader = vtkSmartPointer<vtkPolyDataReader>::New();
	clams_reader->SetFileName("../../../clams/CLaMS_Nabro.vtk");
	clams_reader->Update();
	vtkSmartPointer<vtkPolyData> clams_traj_data = clams_reader->GetOutput();
	auto trajs = clams_traj_data->GetCellData()->GetArray("seed_id");
	int cell_num = clams_traj_data->GetNumberOfCells();
	//clams_traj_data->GetPointData()->SetScalars(trajs);

	std::vector<vtkSmartPointer<vtkPolyData>> tropo_vec;
	std::string tropo_file_prefix = "../../../tropopause/Tropoause-VISContest_";
	std::string file_suffix = ".vtk";
	std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, double>>> txya;

	for (int i = 0; i <= 84; ++i) {
		std::string cur_file = tropo_file_prefix + (i > 9 ? "0" + std::to_string(i) : "00" + std::to_string(i)) + file_suffix;
		vtkSmartPointer<vtkStructuredGridReader> tropo_reader = vtkSmartPointer<vtkStructuredGridReader>::New();
		tropo_reader->SetFileName(cur_file.c_str());
		tropo_reader->Update();
		vtkSmartPointer<vtkStructuredGridGeometryFilter> geometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
		geometryFilter->SetInputConnection(tropo_reader->GetOutputPort());
		geometryFilter->Update();
		vtkSmartPointer<vtkPolyData> tropo_data = geometryFilter->GetOutput();
		auto tropo_1 = tropo_data->GetPointData()->GetArray("trop_1");
		auto time = tropo_data->GetPointData()->GetArray("time");
		int pnum = tropo_data->GetNumberOfPoints();
		//cout << pnum << endl;
		
		for (int i = 0; i < pnum; ++i) {
			int tm = time->GetTuple1(i);
			cout.setf(ios::fixed, ios::floatfield);
			if (i == 0) cout << tm << endl;
			double xyz[3];
			tropo_data->GetPoints()->GetPoint(i, xyz);
			txya[tm][xyz[0]][xyz[1]] = tropo_1->GetTuple1(i);
		}
		tropo_data->GetPointData()->SetScalars(tropo_1);
		tropo_vec.push_back(tropo_data);
	}

	vtkSmartPointer<vtkPoints> stratos_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> stratos_colors = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkCellArray> stratos_cells = vtkSmartPointer<vtkCellArray>::New();
	
	stratos_colors->SetNumberOfValues(clams_traj_data->GetNumberOfPoints());
	for (int i = 0; i < clams_traj_data->GetNumberOfPoints(); ++i) {
		stratos_colors->SetValue(i, 0);
	}
	

	auto altitude = clams_traj_data->GetPointData()->GetArray("altitude");
	auto time = clams_traj_data->GetPointData()->GetArray("time");
	int newid = 0, oldid = 0;

	for (int cn = 0; cn < cell_num; ++cn) {
		auto cur_cell = clams_traj_data->GetCell(cn);
		int seed_id = trajs->GetTuple1(cn);
		//cout << "seed id: " << seed_id << endl;
		bool flag = false;
		int pnum = cur_cell->GetNumberOfPoints();
		for (int pn = 0; pn < pnum; ++pn) {
			int pid = cur_cell->GetPointId(pn);
			double xyz[3];
			clams_traj_data->GetPoints()->GetPoint(pid, xyz);
			double x = xyz[0], y = xyz[1];
			if (abs(x - nabro_x) < 3 && abs(y - nabro_y) < 3) {
				flag = true;
				break;
			}
		}
		if (!flag) {
			for (int pn = 0; pn < pnum; ++pn) {
				int pid = cur_cell->GetPointId(pn);
				stratos_colors->SetValue(pid, 0);
			}
			continue;
		}
		std::vector<std::pair<double, int>> pts;
		double mintm = 1 << 30, minid = 0;
		for (int pn = 0; pn < pnum; ++pn) {
			int pid = cur_cell->GetPointId(pn);
			double alti = altitude->GetTuple1(pid);
			double tm = time->GetTuple1(pid);
			if(tm > 360460800 + 216 * 3600 + 8*24*3600){stratos_colors->SetValue(pid,0);continue;}
			double xyz[3];
			clams_traj_data->GetPoints()->GetPoint(pid, xyz);
			double x = xyz[0], y = xyz[1];
			int tmkey = static_cast<int>(tm) / 100 * 100;
			//cout.setf(ios::fixed, ios::floatfield);
			//if (pn == 0) cout << "cell id: " << cn<< " time: "<<tmkey << endl;
			int xkey = static_cast<int>(x);
			int ykey = static_cast<int>(y);
			if (txya.find(tmkey) == txya.end() ||
				txya[tmkey].find(xkey) == txya[tmkey].end() ||
				txya[tmkey][xkey].find(ykey) == txya[tmkey][xkey].end()) {
				stratos_colors->SetValue(pid, 0);
				continue;
			}
			double stratos_alti = txya[tmkey][xkey][ykey];
			if (alti >= stratos_alti && tm < mintm) {
				mintm = tm;
				minid = pid;
			}
			// if (tm > 360460800 + 216 * 3600) {
			// 	stratos_colors->SetValue(pid, 0);
			// 	continue;
			// }
			// double xyz[3];
			// clams_traj_data->GetPoints()->GetPoint(pid, xyz);
			// double x = xyz[0], y = xyz[1];
			// int tmkey = static_cast<int>(tm) / 100 * 100;
			// //cout.setf(ios::fixed, ios::floatfield);
			// //if (pn == 0) cout << "cell id: " << cn<< " time: "<<tmkey << endl;
			// int xkey = static_cast<int>(x);
			// int ykey = static_cast<int>(y);
			// if (xkey < -160 || xkey > 160) {
			// 	stratos_colors->SetValue(pid, 0);
			// 	continue;
			// }
			// if (txya.find(tmkey) == txya.end() ||
			// 	txya[tmkey].find(xkey) == txya[tmkey].end() ||
			// 	txya[tmkey][xkey].find(ykey) == txya[tmkey][xkey].end()) {
			// 	stratos_colors->SetValue(pid, 0);
			// 	continue;
			// }
			// //cout << "newpoint\n";
			// double stratos_alti = txya[tmkey][xkey][ykey];
			// if (alti < stratos_alti) {
			// 	stratos_colors->SetValue(pid, 0);
			// 	continue;
			// }
			//cout << "stratospoint\n";
			// stratos_colors->SetValue(pid, seed_id);
		}
		stratos_colors->SetValue(minid, seed_id);
	}
	
	clams_traj_data->GetPointData()->SetScalars(stratos_colors);
	

	//vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//mapper->SetInputConnection(clams_reader->GetOutputPort());
	//mapper->ScalarVisibilityOff();
	vtkSmartPointer<vtkLookupTable> tropolut = vtkSmartPointer<vtkLookupTable>::New();
	/*tropolut->SetNumberOfTableValues(256);
	tropolut->SetHueRange(0.0, 0.667);
	tropolut->SetAlphaRange(0.0, 1.0);
	tropolut->SetTableRange(0, 700);*/
	tropolut->SetScaleToLinear();
	tropolut->SetNumberOfTableValues(numColors);
	double r, g, b;
	for (int i = 0; i < numColors; i++)
	{
		double val = min + ((double)i / numColors) * range;
		getColorCorrespondingTovalue(val, r, g, b);
		
		if (r==1 && g==1 && b==1) tropolut->SetTableValue(i, r, g, b, 0.0);
		else 
			tropolut->SetTableValue(i, r, g, b, 1.0);
	}
	tropolut->Build();

	vtkSmartPointer<vtkPolyDataMapper> tropomapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//tropomapper->SetInputConnection(delaunay->GetOutputPort());
	tropomapper->SetInputData(clams_traj_data);
	tropomapper->SetLookupTable(tropolut);
	//tropomapper->UseLookupTableScalarRangeOn();
	tropomapper->SetScalarRange(0, 650);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	//actor->GetProperty()->SetOpacity(1);
	//actor->GetProperty()->SetEdgeVisibility(true);
	actor->SetMapper(tropomapper);

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

int main()
{
	readVTKFile();
	return 0;
}
