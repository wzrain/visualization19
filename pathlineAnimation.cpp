#include <vtkAnimationCue.h>
#include <vtkSphereSource.h>
#include <vtkCommand.h>
#include <vtkAnimationScene.h>
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
#include <unordered_set>

// Nabro:                 13'22'' N    41'42'' E
// Grímsvötn:             64°25′12″N   17°19′48″W
// Puyehue-Cordón Caulle: 40°35′25″S   72°07′02″W

static const double min = 0;
static const double max = 650;
static const double range = max-min;
static const double numColors = 100;
static const long timeStart = 360460800; // P
static const long timeEnd = 367664384; // P
// static const long timeStart = 361217800; // N
// static const long timeEnd = 368582384; // N
static const double X = -72.0f; // P
static const double Y = -40.0f; // P
// static const double X = 41.0f; // N
// static const double Y = 13.0f; // N
static const double delta = 1.0f;
static const double timeDelta = 3600 * 6;
static const double timeGap = 3600 * 24; // one day
static const int animationStartTime = 0;
static const int animationEndTime = 100;
static const int animationTime = animationEndTime - animationStartTime;

struct P{
  double trajId;
  double x, y;
  vtkIdType pointId;
  P(void){}
  P(double ti, double xx, double yy, vtkIdType pi){trajId = ti; x = xx; y = yy; pointId = pi;}
};

class CueAnimator
{
public:
  CueAnimator(
    //vtkPolyDataMapper *m = 0, 
    std::vector<std::vector<P>> pI,
    vtkPolyDataReader *r = 0, 
    vtkActor *a = 0, 
    vtkScalarBarActor *l = 0, 
    vtkCamera *c=0, 
    vtkFloatArray *s = 0)
  {
      this->Reader=r;
      // this->Mapper=m;
      this->Actor=a;
      this->Legend=l;
      this->Camera=c;
      this->Scale=s;
      this->pointInfo=pI;
  }

  ~CueAnimator()
  {
      this->Cleanup();
  }

  void StartCue(vtkAnimationCue::AnimationCueInfo *vtkNotUsed(info),
                vtkRenderer *renderer)
  {
      cout << "*** IN StartCue " << endl;
      // this->SphereSource=vtkSphereSource::New();
      // this->SphereSource->SetRadius(0.5);

      // this->Mapper=vtkPolyDataMapper::New();
      // this->Mapper->SetInputConnection(this->SphereSource->GetOutputPort());

      // this->Actor=vtkActor::New();
      // this->Actor->SetMapper(this->Mapper);

      renderer->AddActor(this->Actor);
      renderer->AddActor2D(this->Legend);
      renderer->SetActiveCamera(this->Camera);
      renderer->ResetCamera();
      renderer->SetBackground(.5, .5, .5);
      renderer->ResetCameraClippingRange();
      renderer->Render();

      // ren->AddActor(this->Actor);
      // ren->ResetCamera();
      // ren->Render();
      cout << "*** OUT StartCue " << endl;
  }

  void Tick(vtkAnimationCue::AnimationCueInfo *info,
            vtkRenderer *render)
  {
      cout << "*** IN Tick " << endl;
      int idx = static_cast<double>(info->AnimationTime - info->StartTime);
      std::cout << idx << ' ' << pointInfo[idx].size() << std::endl;  
      for(int i = 0; i < pointInfo[idx].size(); i++){
        // double px = pointInfo[idx][i].x;
        // double py = pointInfo[idx][i].y;
        // if()
        Scale->SetValue(pointInfo[idx][i].pointId,pointInfo[idx][i].trajId);
      }
      Reader->GetOutput()->GetPointData()->SetScalars(Scale);
      Reader->Update();
      render->Render();
      cout << "*** OUT Tick " << endl;
  }

  void EndCue(vtkAnimationCue::AnimationCueInfo *vtkNotUsed(info),
              vtkRenderer *render)
  {
      cout << "*** IN EndCue " << endl;
      (void)render;
      this->Cleanup();
      cout << "*** OUT EndCue " << endl;
  }

protected:
  vtkPolyDataReader *Reader;
  // vtkPolyDataMapper *Mapper;
  vtkActor *Actor;
  vtkScalarBarActor *Legend;
  vtkCamera *Camera;
  vtkFloatArray *Scale;
  std::vector<std::vector<P>> pointInfo;

  void Cleanup()
  {

      if(this->Reader!=0)
      {
        this->Reader->Delete();
        this->Reader=0;
      }
      if(this->Actor!=0)
      {
        this->Actor->Delete();
        this->Actor=0;
      }
      if(this->Legend!=0)
      {
        this->Legend->Delete();
        this->Legend=0;
      }
      if(this->Camera!=0)
      {
        this->Camera->Delete();
        this->Camera=0;
      }
      if(this->Scale!=0)
      {
        this->Scale->Delete();
        this->Scale=0;
      }
  }
};

class vtkAnimationCueObserver : public vtkCommand
{
public:
  static vtkAnimationCueObserver *New()
  {
      return new vtkAnimationCueObserver;
  }

  virtual void Execute(vtkObject *vtkNotUsed(caller),
                       unsigned long event,
                       void *calldata)
  {
      cout << "*** IN Exe " << endl;
      if(this->Animator!=0 && this->Renderer!=0)
      {
        vtkAnimationCue::AnimationCueInfo *info=
          static_cast<vtkAnimationCue::AnimationCueInfo *>(calldata);
        switch(event)
        {
          case vtkCommand::StartAnimationCueEvent:
            this->Animator->StartCue(info,this->Renderer);
            break;
          case vtkCommand::EndAnimationCueEvent:
            this->Animator->EndCue(info,this->Renderer);
            break;
          case vtkCommand::AnimationCueTickEvent:
            this->Animator->Tick(info,this->Renderer);
            break;
        }
      }
      if(this->RenWin!=0)
      {
        this->RenWin->Render();
      }
  }

  vtkRenderer *Renderer;
  vtkRenderWindow *RenWin;
  CueAnimator *Animator;
protected:
  vtkAnimationCueObserver()
  {
      this->Renderer=0;
      this->Animator=0;
      this->RenWin=0;
  }
};

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

void readVTKFile(void)
{
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renWinInteractor->SetRenderWindow(renWin);
  renWin->AddRenderer(renderer);
  renWin->Render();

  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  // reader->SetFileName("../../../clams/CLaMS_Nabro.vtk");
  reader->SetFileName("../../../clams/CLaMS_Puyehue.vtk");
  reader->Update();

  vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
  std::cout << "Cell Number: " << data->GetNumberOfCells() << std::endl;
  std::cout << "Point Number: " << data->GetNumberOfPoints() << std::endl;

  // done for scale
  vtkSmartPointer<vtkFloatArray> scale = vtkSmartPointer<vtkFloatArray>::New();
  scale->SetNumberOfValues(data->GetNumberOfPoints());
  for(int i = 0; i < data->GetNumberOfPoints(); i++){scale->SetValue(i,0);}

  // save point data
  auto time = data->GetPointData()->GetArray("time");
  auto trajs = data->GetCellData()->GetArray("seed_id"); // the cell type is VTK_POLY_LINE
  long timeThresholdUp = timeStartHour(5*24);
  long timeThresholdDown = timeStartHour(0);
  std::vector<std::vector<P>> pointInfo(animationTime+1);
  // std::unordered_set<vtkIdType> visited;
  for(int i = 0; i < trajs->GetNumberOfTuples(); i++){
    int pointNum = data->GetCell(i)->GetPointIds()->GetNumberOfIds();
    std::vector<double> PTime;
    std::vector<double> TIdx;
    std::vector<vtkIdType> PIdx;
    std::vector<double> PX;
    std::vector<double> PY;
    
    bool near = false;
    for(int j = 0; j < pointNum; j++){
      float pTime = time->GetTuple1(data->GetCell(i)->GetPointId(j)); 
      double co[3];
      data->GetPoints()->GetPoint(data->GetCell(i)->GetPointId(j),co);
      if(pTime < timeStart + 20 * timeGap){
        if(!near and std::abs(co[0]-X) < delta and std::abs(co[1]-Y) < delta 
          and pTime < timeThresholdDown + timeDelta
          ){near = true;}
        PTime.push_back(pTime);
        TIdx.push_back(trajs->GetTuple1(i));
        PIdx.push_back(data->GetCell(i)->GetPointId(j));
        PX.push_back(co[0]);
        PY.push_back(co[1]);
      }
    }
    if(near){
      for(int j = 0; j < PTime.size(); j++){
        if(PX[j] > -170 and PX[j] < 170){
          double ptime = PTime[j];
          int idx = (ptime - timeStart) / timeGap;
          if(idx < pointInfo.size()){
            pointInfo[idx].push_back(P(TIdx[j],PX[j],PY[j],PIdx[j]));
          }
          // idx->SetValue(PIdx[j],PTime[j]);
        }
      }
    }
  }
  // data->GetPointData()->SetScalars(idx);
  
  // lookuptable
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

  // legend done
  vtkSmartPointer <vtkScalarBarActor> legend = vtkScalarBarActor::New();
  legend->SetLookupTable(lookupTable);
  legend->SetNumberOfLabels(5);
  legend->SetTitle("id of trajectory");
  legend->SetVerticalTitleSeparation(6);
  legend->GetPositionCoordinate()->SetValue(0.88,0.1);
  legend->SetWidth(0.05);

  // mapper done
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SetLookupTable(lookupTable);
  mapper->SetScalarRange(0,650);

  // actor done
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // camera done
  vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
  aCamera->SetViewUp(0, 1, 0);
  aCamera->SetPosition(0, 0, 1);
  aCamera->SetFocalPoint(0, 0, 0);
  aCamera->ComputeViewPlaneNormal();

  vtkSmartPointer<vtkAnimationScene> scene = vtkSmartPointer<vtkAnimationScene>::New();
  scene->SetModeToRealTime();
  // scene->SetModeToSequence();
  scene->SetLoop(0);
  scene->SetFrameRate(1);
  scene->SetStartTime(3);
  scene->SetEndTime(80);

  vtkSmartPointer<vtkAnimationCue> cue = vtkSmartPointer<vtkAnimationCue>::New();
  cue->SetStartTime(5);
  cue->SetEndTime(90);
  scene->AddCue(cue);

  CueAnimator animator(pointInfo,reader,actor,legend,aCamera,scale);

  vtkSmartPointer<vtkAnimationCueObserver> observer = vtkSmartPointer<vtkAnimationCueObserver>::New();
  observer->Renderer=renderer;
  observer->Animator=&animator;
  observer->RenWin=renWin;

  cue->AddObserver(vtkCommand::StartAnimationCueEvent,observer);
  cue->AddObserver(vtkCommand::EndAnimationCueEvent,observer);
  cue->AddObserver(vtkCommand::AnimationCueTickEvent,observer);

  scene->Play();
  scene->Stop();

  // renderer->AddActor(actor);
  // renderer->AddActor2D(legend);
  // renderer->SetActiveCamera(aCamera);
  // renderer->ResetCamera();
  // renderer->SetBackground(.5, .5, .5);
  // renderer->ResetCameraClippingRange();

  // renWin->SetSize(600, 600);
  // renWin->Render();
  // renWinInteractor->Initialize();
  renWinInteractor->Start();
}

int main(int argc, char *argv[]){
  readVTKFile();
  return 0;
}



// class CueAnimator
// {
// public:
//   CueAnimator()
//   {
//       this->SphereSource=0;
//       this->Mapper=0;
//       this->Actor=0;
//   }

//   ~CueAnimator()
//   {
//       this->Cleanup();
//   }

//   void StartCue(vtkAnimationCue::AnimationCueInfo *vtkNotUsed(info),
//                 vtkRenderer *ren)
//   {
//       cout << "*** IN StartCue " << endl;
//       this->SphereSource=vtkSphereSource::New();
//       this->SphereSource->SetRadius(0.5);

//       this->Mapper=vtkPolyDataMapper::New();
//       this->Mapper->SetInputConnection(this->SphereSource->GetOutputPort());

//       this->Actor=vtkActor::New();
//       this->Actor->SetMapper(this->Mapper);

//       ren->AddActor(this->Actor);
//       ren->ResetCamera();
//       ren->Render();
//   }

//   void Tick(vtkAnimationCue::AnimationCueInfo *info,
//             vtkRenderer *ren)
//   {
//       cout << "*** IN Tick " << endl;
//       double newradius=0.1 +
//         (static_cast<double>(info->AnimationTime -
//                              info->StartTime)/
//          static_cast<double>(info->EndTime-info->StartTime)) * 1;
//       this->SphereSource->SetRadius(newradius);
//       this->SphereSource->Update();
//       ren->Render();
//   }

//   void EndCue(vtkAnimationCue::AnimationCueInfo *vtkNotUsed(info),
//               vtkRenderer *ren)
//   {
//       cout << "*** IN EndCue " << endl;
//       (void)ren;
//       // don't remove the actor for the regression image.
// //      ren->RemoveActor(this->Actor);
//       this->Cleanup();
//   }

// protected:
//   vtkSphereSource *SphereSource;
//   vtkPolyDataMapper *Mapper;
//   vtkActor *Actor;

//   void Cleanup()
//   {
//       if(this->SphereSource!=0)
//       {
//         this->SphereSource->Delete();
//         this->SphereSource=0;
//       }

//       if(this->Mapper!=0)
//       {
//         this->Mapper->Delete();
//         this->Mapper=0;
//       }
//       if(this->Actor!=0)
//       {
//         this->Actor->Delete();
//         this->Actor=0;
//       }
//   }
// };

// class vtkAnimationCueObserver : public vtkCommand
// {
// public:
//   static vtkAnimationCueObserver *New()
//   {
//       return new vtkAnimationCueObserver;
//   }

//   virtual void Execute(vtkObject *vtkNotUsed(caller),
//                        unsigned long event,
//                        void *calldata)
//   {
//       cout << "*** IN Exe " << endl;
//       if(this->Animator!=0 && this->Renderer!=0)
//       {
//         vtkAnimationCue::AnimationCueInfo *info=
//           static_cast<vtkAnimationCue::AnimationCueInfo *>(calldata);
//         switch(event)
//         {
//           case vtkCommand::StartAnimationCueEvent:
//             this->Animator->StartCue(info,this->Renderer);
//             break;
//           case vtkCommand::EndAnimationCueEvent:
//             this->Animator->EndCue(info,this->Renderer);
//             break;
//           case vtkCommand::AnimationCueTickEvent:
//             this->Animator->Tick(info,this->Renderer);
//             break;
//         }
//       }
//       if(this->RenWin!=0)
//       {
//         this->RenWin->Render();
//       }
//   }

//   vtkRenderer *Renderer;
//   vtkRenderWindow *RenWin;
//   CueAnimator *Animator;
// protected:
//   vtkAnimationCueObserver()
//   {
//       this->Renderer=0;
//       this->Animator=0;
//       this->RenWin=0;
//   }
// };

// int main(int, char *[])
// {
//   // Create the graphics structure. The renderer renders into the
//   // render window.
//   vtkSmartPointer<vtkRenderWindowInteractor> iren =
//     vtkSmartPointer<vtkRenderWindowInteractor>::New();
//   vtkSmartPointer<vtkRenderer> ren1 =
//     vtkSmartPointer<vtkRenderer>::New();
//   vtkSmartPointer<vtkRenderWindow> renWin =
//     vtkSmartPointer<vtkRenderWindow>::New();
//   renWin->SetMultiSamples(0);
//   iren->SetRenderWindow(renWin);
//   renWin->AddRenderer(ren1);
//   renWin->Render();

//   // Create an Animation Scene
//   vtkSmartPointer<vtkAnimationScene> scene =
//     vtkSmartPointer<vtkAnimationScene>::New();

//   scene->SetModeToRealTime();
//   //scene->SetModeToSequence();

//   scene->SetLoop(0);
//   scene->SetFrameRate(5);
//   scene->SetStartTime(3);
//   scene->SetEndTime(20);

//   // Create an Animation Cue.
//   vtkSmartPointer<vtkAnimationCue> cue1 =
//     vtkSmartPointer<vtkAnimationCue>::New();
//   cue1->SetStartTime(5);
//   cue1->SetEndTime(23);
//   scene->AddCue(cue1);

//   // Create cue animator;
//   CueAnimator animator;

//   // Create Cue observer.
//   vtkSmartPointer<vtkAnimationCueObserver> observer =
//     vtkSmartPointer<vtkAnimationCueObserver>::New();
//   observer->Renderer=ren1;
//   observer->Animator=&animator;
//   observer->RenWin=renWin;

//   cue1->AddObserver(vtkCommand::StartAnimationCueEvent,observer);
//   cue1->AddObserver(vtkCommand::EndAnimationCueEvent,observer);
//   cue1->AddObserver(vtkCommand::AnimationCueTickEvent,observer);

//   scene->Play();
//   scene->Stop();

//   iren->Start();

//   return EXIT_SUCCESS;
// }
