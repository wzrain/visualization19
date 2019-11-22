## airs
### Reader example
```c++
vtkSmartPointer<vtkPolyDataReader> reader =vtkSmartPointer<vtkPolyDataReader>::New();
reader->SetFileName("...\\airs\\volcano_2011_150_am.vtk");
reader->Update();

vtkSmartPointer<vtkPolyData> data = reader->GetOutput();
auto ashIndex = data->GetPointData()->GetArray("SO2");
data->GetPointData()->SetScalars(ashIndex);
```
### Attribute name
 - time
 - ash
 - SO2

## airs_resampled
### Reader example
```c++
```
### Attribute name

## clams
### Reader example
Same as airs

### Attribute name
 - time
 - altitude
 - pressure
 - temperature
 - pot_temperature
 - pot_vorticity

## mipas
### Reader example
Same as airs

### Attribute name
 - time
 - altitude
 - orbit_id
 - profile_id
 - detection