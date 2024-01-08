#include "CompuFluidDyna.hpp"


// Standard lib
#include <numbers>
#include <vector>

// GLUT lib
#include "freeglut/include/GL/freeglut.h"

// Algo headers
#include "FileIO/FileInput.hpp"
#include "Math/Field.hpp"
#include "Math/Vec.hpp"
#include "Util/Random.hpp"

// Global headers
#include "Data.hpp"


// Link to shared sandbox data
extern Data D;


void CompuFluidDyna::SetUpUIData() {
  // Draw the scatter data
  const int yCursor= std::min(std::max((int)std::round((float)(nY - 1) * D.UI[SlicePlotY______].F()), 0), nY - 1);
  const int zCursor= std::min(std::max((int)std::round((float)(nZ - 1) * D.UI[SlicePlotZ______].F()), 0), nZ - 1);
  D.scatLegend.resize(4);
  D.scatLegend[0]= "Horiz VZ";
  D.scatLegend[1]= "Verti VY";
  D.scatLegend[2]= "Horiz P";
  D.scatLegend[3]= "Verti P";
  D.scatData.resize(4);
  for (int k= 0; k < (int)D.scatData.size(); k++)
    D.scatData[k].clear();
  if (nZ > 1) {
    for (int y= 0; y < nY; y++) {
      D.scatData[0].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), VelZ[nX / 2][y][zCursor]}));
      D.scatData[2].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), Pres[nX / 2][y][zCursor]}));
    }
  }
  if (nY > 1) {
    for (int z= 0; z < nZ; z++) {
      D.scatData[1].push_back(std::array<double, 2>({VelY[nX / 2][yCursor][z], (double)z / (double)(nZ - 1)}));
      D.scatData[3].push_back(std::array<double, 2>({Pres[nX / 2][yCursor][z], (double)z / (double)(nZ - 1)}));
    }
  }

  // Add hard coded experimental values for lid driven cavity flow benchmark
  if (D.UI[Scenario________].I() == 3) {
    // Clear unnecessary scatter data
    D.scatData[2].clear();
    D.scatData[3].clear();
    // Allocate required scatter data
    D.scatLegend.resize(8);
    D.scatLegend[4]= "Ghia Re1k";
    D.scatLegend[5]= "Ghia Re1k";
    D.scatLegend[6]= "Ertu Re1k";
    D.scatLegend[7]= "Ertu Re1k";
    D.scatData.resize(8);
    D.scatData[4].clear();
    D.scatData[5].clear();
    D.scatData[6].clear();
    D.scatData[7].clear();
    // TODO add vorticity data from https://www.acenumerics.com/the-benchmarks.html
    // Data from Ghia 1982 http://www.msaidi.ir/upload/Ghia1982.pdf
    const std::vector<double> GhiaData0X({0.0000, +0.0625, +0.0703, +0.0781, +0.0938, +0.1563, +0.2266, +0.2344, +0.5000, +0.8047, +0.8594, +0.9063, +0.9453, +0.9531, +0.9609, +0.9688, +1.0000});  // coord on horiz slice
    const std::vector<double> GhiaData1Y({0.0000, +0.0547, +0.0625, +0.0703, +0.1016, +0.1719, +0.2813, +0.4531, +0.5000, +0.6172, +0.7344, +0.8516, +0.9531, +0.9609, +0.9688, +0.9766, +1.0000});  // coord on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.0923, +0.1009, +0.1089, +0.1232, +0.1608, +0.1751, +0.1753, +0.0545, -0.2453, -0.2245, -0.1691, -0.1031, -0.0886, -0.0739, -0.0591, +0.0000});  // Re 100   verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.0372, -0.0419, -0.0478, -0.0643, -0.1015, -0.1566, -0.2109, -0.2058, -0.1364, +0.0033, +0.2315, +0.6872, +0.7372, +0.7887, +0.8412, +1.0000});  // Re 100   horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.1836, +0.1971, +0.2092, +0.2297, +0.2812, +0.3020, +0.3017, +0.0519, -0.3860, -0.4499, -0.2383, -0.2285, -0.1925, -0.1566, -0.1215, +0.0000});  // Re 400   verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.0819, -0.0927, -0.1034, -0.1461, -0.2430, -0.3273, -0.1712, -0.1148, +0.0214, +0.1626, +0.2909, +0.5589, +0.6176, +0.6844, +0.7584, +1.0000});  // Re 400   horiz vel on verti slice
    const std::vector<double> GhiaData0Y({0.0000, +0.2749, +0.2901, +0.3035, +0.3263, +0.3710, +0.3308, +0.3224, +0.0253, -0.3197, -0.4267, -0.5150, -0.3919, -0.3371, -0.2767, -0.2139, +0.0000});  // Re 1000  verti vel on horiz slice
    const std::vector<double> GhiaData1X({0.0000, -0.1811, -0.2020, -0.2222, -0.2973, -0.3829, -0.2781, -0.1065, -0.0608, +0.0570, +0.1872, +0.3330, +0.4660, +0.5112, +0.5749, +0.6593, +1.0000});  // Re 1000  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.3956, +0.4092, +0.4191, +0.4277, +0.3712, +0.2903, +0.2819, +0.0100, -0.3118, -0.3740, -0.4431, -0.5405, -0.5236, -0.4743, -0.3902, +0.0000});  // Re 3200  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.3241, -0.3534, -0.3783, -0.4193, -0.3432, -0.2443, -0.8664, -0.0427, +0.0716, +0.1979, +0.3468, +0.4610, +0.4655, +0.4830, +0.5324, +1.0000});  // Re 3200  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4245, +0.4333, +0.4365, +0.4295, +0.3537, +0.2807, +0.2728, +0.0095, -0.3002, -0.3621, -0.4144, -0.5288, -0.5541, -0.5507, -0.4977, +0.0000});  // Re 5000  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4117, -0.4290, -0.4364, -0.4044, -0.3305, -0.2286, -0.0740, -0.0304, +0.0818, +0.2009, +0.3356, +0.4604, +0.4599, +0.4612, +0.4822, +1.0000});  // Re 5000  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4398, +0.4403, +0.4356, +0.4182, +0.3506, +0.2812, +0.2735, +0.0082, -0.3045, -0.3621, -0.4105, -0.4859, -0.5235, -0.5522, -0.5386, +0.0000});  // Re 7500  verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4315, -0.4359, -0.4303, -0.3832, -0.3239, -0.2318, -0.0750, -0.0380, +0.0834, +0.2059, +0.3423, +0.4717, +0.4732, +0.4705, +0.4724, +1.0000});  // Re 7500  horiz vel on verti slice
    // const std::vector<double> GhiaData0Y({0.0000, +0.4398, +0.4373, +0.4312, +0.4149, +0.3507, +0.2800, +0.2722, +0.0083, -0.3072, -0.3674, -0.4150, -0.4586, -0.4910, -0.5299, -0.5430, +0.0000});  // Re 10000 verti vel on horiz slice
    // const std::vector<double> GhiaData1X({0.0000, -0.4274, -0.4254, -0.4166, -0.3800, -0.3271, -0.2319, -0.0754, +0.0311, +0.0834, +0.2067, +0.3464, +0.4780, +0.4807, +0.4778, +0.4722, +1.0000});  // Re 10000 horiz vel on verti slice
    // Data from Erturk 2005 https://arxiv.org/pdf/physics/0505121.pdf
    const std::vector<double> ErtuData0X({0.0000, +0.0150, +0.0300, +0.0450, +0.0600, +0.0750, +0.0900, +0.1050, +0.1200, +0.1350, +0.1500, +0.5000, +0.8500, +0.8650, +0.8800, +0.8950, +0.9100, +0.9250, +0.9400, +0.9550, +0.9700, +0.9850, +1.0000});  // coord on horiz slice
    const std::vector<double> ErtuData1Y({0.0000, +0.0200, +0.0400, +0.0600, +0.0800, +0.1000, +0.1200, +0.1400, +0.1600, +0.1800, +0.2000, +0.5000, +0.9000, +0.9100, +0.9200, +0.9300, +0.9400, +0.9500, +0.9600, +0.9700, +0.9800, +0.9900, +1.0000});  // coord on verti slice
    const std::vector<double> ErtuData0Y({0.0000, +0.1019, +0.1792, +0.2349, +0.2746, +0.3041, +0.3273, +0.3460, +0.3605, +0.3705, +0.3756, +0.0258, -0.4028, -0.4407, -0.4803, -0.5132, -0.5263, -0.5052, -0.4417, -0.3400, -0.2173, -0.0973, +0.0000});  // Re 1000  verti vel on horiz slice
    const std::vector<double> ErtuData1X({0.0000, -0.0757, -0.1392, -0.1951, -0.2472, -0.2960, -0.3381, -0.3690, -0.3854, -0.3869, -0.3756, -0.0620, +0.3838, +0.3913, +0.3993, +0.4101, +0.4276, +0.4582, +0.5102, +0.5917, +0.7065, +0.8486, +1.0000});  // Re 1000  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.1607, +0.2633, +0.3238, +0.3649, +0.3950, +0.4142, +0.4217, +0.4187, +0.4078, +0.3918, +0.0160, -0.3671, -0.3843, -0.4042, -0.4321, -0.4741, -0.5268, -0.5603, -0.5192, -0.3725, -0.1675, +0.0000});  // Re 2500  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.1517, -0.2547, -0.3372, -0.3979, -0.4250, -0.4200, -0.3965, -0.3688, -0.3439, -0.3228, -0.0403, +0.4141, +0.4256, +0.4353, +0.4424, +0.4470, +0.4506, +0.4607, +0.4971, +0.5924, +0.7704, +1.0000});  // Re 2500  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2160, +0.3263, +0.3868, +0.4258, +0.4426, +0.4403, +0.4260, +0.4070, +0.3878, +0.3699, +0.0117, -0.3624, -0.3806, -0.3982, -0.4147, -0.4318, -0.4595, -0.5139, -0.5700, -0.5019, -0.2441, +0.0000});  // Re 5000  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2223, -0.3480, -0.4272, -0.4419, -0.4168, -0.3876, -0.3652, -0.3467, -0.3285, -0.3100, -0.0319, +0.4155, +0.4307, +0.4452, +0.4582, +0.4683, +0.4738, +0.4739, +0.4749, +0.5159, +0.6866, +1.0000});  // Re 5000  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2509, +0.3608, +0.4210, +0.4494, +0.4495, +0.4337, +0.4137, +0.3950, +0.3779, +0.3616, +0.0099, -0.3574, -0.3755, -0.3938, -0.4118, -0.4283, -0.4443, -0.4748, -0.5434, -0.5550, -0.2991, +0.0000});  // Re 7500  verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2633, -0.3980, -0.4491, -0.4284, -0.3978, -0.3766, -0.3587, -0.3406, -0.3222, -0.3038, -0.0287, +0.4123, +0.4275, +0.4431, +0.4585, +0.4723, +0.4824, +0.4860, +0.4817, +0.4907, +0.6300, +1.0000});  // Re 7500  horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2756, +0.3844, +0.4409, +0.4566, +0.4449, +0.4247, +0.4056, +0.3885, +0.3722, +0.3562, +0.0088, -0.3538, -0.3715, -0.3895, -0.4078, -0.4256, -0.4411, -0.4592, -0.5124, -0.5712, -0.3419, +0.0000});  // Re 10000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.2907, -0.4259, -0.4469, -0.4142, -0.3899, -0.3721, -0.3543, -0.3361, -0.3179, -0.2998, -0.0268, +0.4095, +0.4243, +0.4398, +0.4556, +0.4711, +0.4843, +0.4917, +0.4891, +0.4837, +0.5891, +1.0000});  // Re 10000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.2940, +0.4018, +0.4522, +0.4563, +0.4383, +0.4180, +0.4004, +0.3840, +0.3678, +0.3519, +0.0080, -0.3508, -0.3682, -0.3859, -0.4040, -0.4221, -0.4388, -0.4534, -0.4899, -0.5694, -0.3762, +0.0000});  // Re 12500 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3113, -0.4407, -0.4380, -0.4054, -0.3859, -0.3685, -0.3506, -0.3326, -0.3146, -0.2967, -0.0256, +0.4070, +0.4216, +0.4366, +0.4523, +0.4684, +0.4833, +0.4937, +0.4941, +0.4833, +0.5587, +1.0000});  // Re 12500 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3083, +0.4152, +0.4580, +0.4529, +0.4323, +0.4132, +0.3964, +0.3801, +0.3641, +0.3483, +0.0074, -0.3481, -0.3654, -0.3828, -0.4005, -0.4186, -0.4361, -0.4505, -0.4754, -0.5593, -0.4041, +0.0000});  // Re 15000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3278, -0.4474, -0.4286, -0.4001, -0.3827, -0.3652, -0.3474, -0.3297, -0.3119, -0.2942, -0.0247, +0.4047, +0.4190, +0.4338, +0.4492, +0.4653, +0.4811, +0.4937, +0.4969, +0.4850, +0.5358, +1.0000});  // Re 15000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3197, +0.4254, +0.4602, +0.4484, +0.4273, +0.4093, +0.3929, +0.3767, +0.3608, +0.3452, +0.0069, -0.3457, -0.3627, -0.3800, -0.3975, -0.4153, -0.4331, -0.4482, -0.4664, -0.5460, -0.4269, +0.0000});  // Re 17500 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3412, -0.4490, -0.4206, -0.3965, -0.3797, -0.3622, -0.3446, -0.3271, -0.3096, -0.2920, -0.0240, +0.4024, +0.4166, +0.4312, +0.4463, +0.4622, +0.4784, +0.4925, +0.4982, +0.4871, +0.5183, +1.0000});  // Re 17500 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3290, +0.4332, +0.4601, +0.4438, +0.4232, +0.4060, +0.3897, +0.3736, +0.3579, +0.3423, +0.0065, -0.3434, -0.3603, -0.3774, -0.3946, -0.4122, -0.4300, -0.4459, -0.4605, -0.5321, -0.4457, +0.0000});  // Re 20000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3523, -0.4475, -0.4143, -0.3936, -0.3769, -0.3595, -0.3422, -0.3248, -0.3074, -0.2899, -0.0234, +0.4001, +0.4142, +0.4287, +0.4436, +0.4592, +0.4754, +0.4906, +0.4985, +0.4889, +0.5048, +1.0000});  // Re 20000 horiz vel on verti slice
    // const std::vector<double> ErtuData0Y({0.0000, +0.3323, +0.4357, +0.4596, +0.4420, +0.4218, +0.4048, +0.3885, +0.3725, +0.3567, +0.3413, +0.0063, -0.3425, -0.3593, -0.3764, -0.3936, -0.4110, -0.4287, -0.4449, -0.4588, -0.5266, -0.4522, +0.0000});  // Re 21000 verti vel on horiz slice
    // const std::vector<double> ErtuData1X({0.0000, -0.3562, -0.4463, -0.4121, -0.3925, -0.3758, -0.3585, -0.3412, -0.3239, -0.3066, -0.2892, -0.0232, +0.3992, +0.4132, +0.4277, +0.4425, +0.4580, +0.4742, +0.4897, +0.4983, +0.4895, +0.5003, +1.0000});  // Re 21000 horiz vel on verti slice
    // Add hard coded experimental values in the scatter plot
    for (int k= 0; k < (int)GhiaData0X.size(); k++) {
      D.scatData[4].push_back(std::array<double, 2>({GhiaData0X[k], GhiaData0Y[k]}));
      D.scatData[5].push_back(std::array<double, 2>({GhiaData1X[k], GhiaData1Y[k]}));
    }
    for (int k= 0; k < (int)ErtuData0X.size(); k++) {
      D.scatData[6].push_back(std::array<double, 2>({ErtuData0X[k], ErtuData0Y[k]}));
      D.scatData[7].push_back(std::array<double, 2>({ErtuData1X[k], ErtuData1Y[k]}));
    }
  }

  // Add hard coded analytical values for Poiseuille flow benchmark
  if (D.UI[Scenario________].I() == 6) {
    // Clear unnecessary scatter data
    D.scatData[0].clear();
    D.scatData[3].clear();
    // Allocate required scatter data
    D.scatLegend.resize(6);
    D.scatLegend[4]= "Analy VY";
    D.scatLegend[5]= "Analy P";
    D.scatData.resize(6);
    D.scatData[4].clear();
    D.scatData[5].clear();
    // Add analytical values in the scatter plot
    const float press0= D.UI[BCPres__________].F();
    const float press1= -D.UI[BCPres__________].F();
    const float kinVisco= D.UI[CoeffDiffuV_____].F();
    if (nY > 1) {
      for (int z= 0; z < nZ; z++) {
        const float width= voxSize * (float)(nY - 1);
        const float height= voxSize * (float)(nZ - 1);
        const float posZ= (float)z * voxSize;
        const float pressDiff= (press1 - press0) / width;
        const float analyVelY= -pressDiff * (1.0f / (2.0f * kinVisco)) * posZ * (height - posZ);
        D.scatData[4].push_back(std::array<double, 2>({analyVelY, (double)z / (double)(nZ - 1)}));
      }
    }
    if (nZ > 1) {
      for (int y= 0; y < nY; y++) {
        const float analyP= press0 + (press1 - press0) * (float)y / (float)(nY - 1);
        D.scatData[5].push_back(std::array<double, 2>({(double)y / (double)(nY - 1), analyP}));
      }
    }
  }

  if (!D.UI[VerboseSolv_____].B()) {
    // Draw the plot data
    D.plotData.resize(5);
    D.plotLegend.resize(5);
    D.plotLegend[0]= "VelMag";
    D.plotLegend[1]= "Smok";
    D.plotLegend[2]= "Pres";
    D.plotLegend[3]= "DiveAbs";
    D.plotLegend[4]= "Vorti";
    if (D.plotData[0].size() < 1000) {
      for (int k= 0; k < (int)D.plotLegend.size(); k++)
        D.plotData[k].push_back(0.0f);
      for (int x= 0; x < nX; x++) {
        for (int y= 0; y < nY; y++) {
          for (int z= 0; z < nZ; z++) {
            D.plotData[0][D.plotData[0].size() - 1]+= std::sqrt(VelX[x][y][z] * VelX[x][y][z] + VelY[x][y][z] * VelY[x][y][z] + VelZ[x][y][z] * VelZ[x][y][z]);
            D.plotData[1][D.plotData[1].size() - 1]+= Smok[x][y][z];
            D.plotData[2][D.plotData[2].size() - 1]+= Pres[x][y][z];
            D.plotData[3][D.plotData[3].size() - 1]+= std::abs(Dive[x][y][z]);
            D.plotData[4][D.plotData[4].size() - 1]+= Vort[x][y][z];
          }
        }
      }
    }
  }
}


void CompuFluidDyna::InitializeScenario() {
  // Get scenario ID and optionnally load bitmap file
  const int scenarioType= D.UI[Scenario________].I();
  const int inputFile= D.UI[InputFile_______].I();
  std::vector<std::vector<std::array<float, 4>>> imageRGBA;
  if (scenarioType == 0) {
    if (inputFile == 0) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/NACA.bmp", imageRGBA, false);
    if (inputFile == 1) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/Nozzle.bmp", imageRGBA, false);
    if (inputFile == 2) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/Pipe.bmp", imageRGBA, false);
    if (inputFile == 3) FileInput::LoadImageBMPFile("./FileInput/FluidScenarios/TeslaValve.bmp", imageRGBA, false);
  }

  // Set scenario values
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Scenario from loaded BMP file
        if (scenarioType == 0 && !imageRGBA.empty()) {
          const float posW= (float)(imageRGBA.size() - 1) * ((float)y + 0.5f) / (float)nY;
          const float posH= (float)(imageRGBA[0].size() - 1) * ((float)z + 0.5f) / (float)nZ;
          const int idxPixelW= std::min(std::max((int)std::round(posW), 0), (int)imageRGBA.size() - 1);
          const int idxPixelH= std::min(std::max((int)std::round(posH), 0), (int)imageRGBA[0].size() - 1);
          const std::array<float, 4> colRGBA= imageRGBA[idxPixelW][idxPixelH];
          if (colRGBA[3] < 0.1f) {
            Solid[x][y][z]= true;
          }
          else {
            if (std::abs(colRGBA[0] - 0.5f) > 0.1f) PreBC[x][y][z]= true;
            if (std::abs(colRGBA[1] - 0.5f) > 0.1f) VelBC[x][y][z]= true;
            if (std::abs(colRGBA[2] - 0.5f) > 0.1f) SmoBC[x][y][z]= true;
          }
          if (PreBC[x][y][z]) {
            PresForced[x][y][z]= D.UI[BCPres__________].F() * ((colRGBA[0] > 0.5f) ? (1.0f) : (-1.0f));
          }
          if (VelBC[x][y][z]) {
            VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
            VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((colRGBA[1] > 0.5f) ? (1.0f) : (-1.0f));
          }
          if (SmoBC[x][y][z]) {
            SmokForced[x][y][z]= D.UI[BCSmok__________].F() * ((colRGBA[2] > 0.5f) ? (1.0f) : (-1.0f));
          }
        }
        // Double opposing inlets
        // |-----------|
        // |           |
        // | (>)   (<) |
        // |           |
        // |-----------|
        if (scenarioType == 1) {
          if (nY > 1 && (y == 0 || y == nY - 1)) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
          else if ((nX > 1 && (x == 0 || x == nX - 1)) ||
                   (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else {
            const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
            const Vec::Vec3<float> posInlet(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
            for (int k= 0; k < 2; k++) {
              if (k == 0 && (posVox - posInlet).norm() > D.UI[ObjectSize0_____].F()) continue;
              if (k == 1 && (posVox - Vec::Vec3<float>(1.0f, 1.0f, 1.0f) + posInlet).norm() > D.UI[ObjectSize1_____].F()) continue;
              VelBC[x][y][z]= true;
              SmoBC[x][y][z]= true;
              VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
              SmokForced[x][y][z]= D.UI[BCSmok__________].F() * ((k == 0) ? (1.0f) : (-1.0f));
            }
          }
        }
        // Circular obstacle in corridor showing vortex shedding
        // Test calib flow separation past cylinder https://link.springer.com/article/10.1007/s00521-020-05079-z
        // ---------------
        // >   O         >
        // ---------------
        if (scenarioType == 2) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
          }
          else if (y == nY - 1) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
          else if (y == 0) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
          else {
            const Vec::Vec3<float> posCell(((float)x + 0.5f) / (float)nX, ((float)y + 0.5f) / (float)nY, ((float)z + 0.5f) / (float)nZ);
            const Vec::Vec3<float> posObstacle(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
            Vec::Vec3<float> dist= (posCell - posObstacle);
            dist[0]*= (float)(nX - 1) * voxSize;
            dist[1]*= (float)(nY - 1) * voxSize;
            dist[2]*= (float)(nZ - 1) * voxSize;
            if (dist.norm() <= std::max(D.UI[ObjectSize0_____].F(), 0.0f))
              Solid[x][y][z]= true;
          }
        }
        // Cavity lid shear benchmark
        //  >>>>>>>>>>>
        // |           |
        // |           |
        // |           |
        // |-----------|
        if (scenarioType == 3) {
          if (y == 0 || y == nY - 1 || z == 0) {
            Solid[x][y][z]= true;
          }
          else if (z == nZ - 1) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
          }
          else if (y == nY / 2 && z > nZ / 2) {
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
        }
        // Flow constriction test with circular hole in a wall
        // -----|||--------
        // >    |||       >
        // >              >
        // >    |||       >
        // -----|||--------
        if (scenarioType == 4) {
          const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
          const Vec::Vec3<float> posWall(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
          const float radHole= D.UI[ObjectSize0_____].F();
          const float radWall= D.UI[ObjectSize1_____].F();
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else if ((posVox - posWall).coeffMul(Vec::Vec3<float>(0.0f, 1.0f, 0.0f)).norm() <= radWall &&
                   (posVox - posWall).coeffMul(Vec::Vec3<float>(1.0f, 0.0f, 1.0f)).norm() > radHole) {
            Solid[x][y][z]= true;
          }
          else if (y == 0) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F();
            VelYForced[x][y][z]= D.UI[BCVelY__________].F();
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
          else if (y == nY - 1) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= 0.0f;
          }
        }
        // Central bloc with initial velocity
        // --------------
        // |    >>>>    |
        // |    >>>>    |
        // --------------
        if (scenarioType == 5) {
          if (((nX == 1) != (std::min(x, nX - 1 - x) > nX / 3)) &&
              ((nY == 1) != (std::min(y, nY - 1 - y) > nY / 3)) &&
              ((nZ == 1) != (std::min(z, nZ - 1 - z) > nZ / 3))) {
            VelX[x][y][z]= D.UI[BCVelX__________].F();
            VelY[x][y][z]= D.UI[BCVelY__________].F();
            VelZ[x][y][z]= D.UI[BCVelZ__________].F();
            Smok[x][y][z]= D.UI[BCSmok__________].F() * ((std::min(z, nZ - 1 - z) < 4 * (nZ - 1) / 9) ? (1.0f) : (-1.0f));
          }
        }
        // Poiseuille/Couette flow in tube with pressure gradient
        // ---------------
        // high        low
        // ---------------
        if (scenarioType == 6) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            VelBC[x][y][z]= true;
            VelXForced[x][y][z]= D.UI[BCVelX__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
            VelYForced[x][y][z]= D.UI[BCVelY__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
            VelZForced[x][y][z]= D.UI[BCVelZ__________].F() * ((z < nZ / 2) ? (-1.0f) : (1.0f));
          }
          else if (nY > 1 && (y == 0 || y == nY - 1)) {
            PreBC[x][y][z]= true;
            PresForced[x][y][z]= D.UI[BCPres__________].F() * ((y < nY / 2) ? (1.0f) : (-1.0f));
          }
          else if (std::max(y, nY - 1 - y) == nY / 2) {
            SmoBC[x][y][z]= true;
            SmokForced[x][y][z]= D.UI[BCSmok__________].F();
          }
        }
        // Thermal convection cell
        // |---cold---|
        // |          |
        // |---warm---|
        if (scenarioType == 7) {
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && (y == 0 || y == nY - 1)) ||
              (nZ > 1 && (z == 0 || z == nZ - 1))) {
            Solid[x][y][z]= true;
          }
          else if (nZ > 1 && std::min(z, nZ - 1 - z) < nZ / 3) {
            Smok[x][y][z]= D.UI[BCSmok__________].F() * ((z > nZ / 2) ? (-1.0f) : (1.0f));
            Smok[x][y][z]+= Random::Val(-0.01f, 0.01f);
          }
        }
        // Pipe of constant diameter with right angle turn and UI parameters to tweak position, curvature, diameters
        // ----------_
        // >          -
        // -----_      |
        //       |     |
        //       |  v  |
        if (scenarioType == 8) {
          const Vec::Vec3<float> posVox= Vec::Vec3<float>(D.boxMin[0], D.boxMin[1], D.boxMin[2]) + voxSize * Vec::Vec3<float>(x + 0.5f, y + 0.5f, z + 0.5f);
          const Vec::Vec3<float> posBend(D.UI[ObjectPosX______].F(), D.UI[ObjectPosY______].F(), D.UI[ObjectPosZ______].F());
          const float radPipe= D.UI[ObjectSize0_____].F();
          const float radBend= D.UI[ObjectSize1_____].F();
          if ((nX > 1 && (x == 0 || x == nX - 1)) ||
              (nY > 1 && y == nY - 1) ||
              (nZ > 1 && z == nZ - 1)) {
            Solid[x][y][z]= true;
          }
          else if (posVox[1] < posBend[1]) {
            if ((posVox - posBend - Vec::Vec3<float>(0.0f, 0.0f, radBend + radPipe)).coeffMul(Vec::Vec3<float>(1.0f, 0.0f, 1.0f)).norm() > radPipe) {
              Solid[x][y][z]= true;
            }
            else if (y == 0) {
              VelBC[x][y][z]= true;
              VelXForced[x][y][z]= D.UI[BCVelX__________].F();
              VelYForced[x][y][z]= D.UI[BCVelY__________].F();
              VelZForced[x][y][z]= D.UI[BCVelZ__________].F();
              // PreBC[x][y][z]= true;
              // PresForced[x][y][z]= D.UI[BCPres__________].F();
              SmoBC[x][y][z]= true;
              SmokForced[x][y][z]= D.UI[BCSmok__________].F();
            }
          }
          else if (posVox[2] < posBend[2]) {
            if ((posVox - posBend - Vec::Vec3<float>(0.0f, radBend + radPipe, 0.0f)).coeffMul(Vec::Vec3<float>(1.0f, 1.0f, 0.0f)).norm() > radPipe) {
              Solid[x][y][z]= true;
            }
            else if (z == 0) {
              PreBC[x][y][z]= true;
              PresForced[x][y][z]= 0.0;
            }
          }
          else if ((posBend + ((posVox - posBend).coeffMul(Vec::Vec3<float>(0.0f, 1.0f, 1.0f)).normalized() * (radBend + radPipe)) - posVox).norm() > radPipe) {
            Solid[x][y][z]= true;
          }
        }
      }
    }
  }
}


// Apply boundary conditions enforcing fixed values to fields
void CompuFluidDyna::ApplyBC(const int iFieldID, std::vector<std::vector<std::vector<float>>>& ioField) {
  // Sweep through the field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Set forced value
        if (Solid[x][y][z] && iFieldID == FieldID::IDSmok) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelX) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelY) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDVelZ) ioField[x][y][z]= 0.0f;
        if (Solid[x][y][z] && iFieldID == FieldID::IDPres) ioField[x][y][z]= 0.0f;
        if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) ioField[x][y][z]= SmokForced[x][y][z] * std::cos(simTime * 2.0f * std::numbers::pi / D.UI[BCSmokTime______].F());
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelX) ioField[x][y][z]= VelXForced[x][y][z];
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelY) ioField[x][y][z]= VelYForced[x][y][z];
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelZ) ioField[x][y][z]= VelZForced[x][y][z];
        if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) ioField[x][y][z]= PresForced[x][y][z];
      }
    }
  }
}


// Addition of one field to an other
void CompuFluidDyna::ImplicitFieldAdd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                      const std::vector<std::vector<std::vector<float>>>& iFieldB,
                                      std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iFieldA[x][y][z] + iFieldB[x][y][z];
}


// Subtraction of one field to an other
void CompuFluidDyna::ImplicitFieldSub(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                      const std::vector<std::vector<std::vector<float>>>& iFieldB,
                                      std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iFieldA[x][y][z] - iFieldB[x][y][z];
}


// Multiplication of field by scalar
void CompuFluidDyna::ImplicitFieldScale(const float iVal,
                                        const std::vector<std::vector<std::vector<float>>>& iField,
                                        std::vector<std::vector<std::vector<float>>>& oField) {
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        oField[x][y][z]= iField[x][y][z] * iVal;
}


// Dot product between two fields
float CompuFluidDyna::ImplicitFieldDotProd(const std::vector<std::vector<std::vector<float>>>& iFieldA,
                                           const std::vector<std::vector<std::vector<float>>>& iFieldB) {
  float val= 0.0f;
  for (int x= 0; x < nX; x++)
    for (int y= 0; y < nY; y++)
      for (int z= 0; z < nZ; z++)
        val+= iFieldA[x][y][z] * iFieldB[x][y][z];
  return val;
}


// Perform a matrix-vector multiplication without explicitly assembling the Laplacian matrix
void CompuFluidDyna::ImplicitFieldLaplacianMatMult(const int iFieldID, const float iTimeStep,
                                                   const bool iDiffuMode, const float iDiffuCoeff, const bool iPrecondMode,
                                                   const std::vector<std::vector<std::vector<float>>>& iField,
                                                   std::vector<std::vector<std::vector<float>>>& oField) {
  // Precompute value
  const float diffuVal= iDiffuCoeff * iTimeStep / (voxSize * voxSize);
  // Sweep through the field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        // Skip solid or fixed values
        if (Solid[x][y][z]) continue;
        if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
        if (VelBC[x][y][z] && (iFieldID == FieldID::IDVelX || iFieldID == FieldID::IDVelY || iFieldID == FieldID::IDVelZ)) continue;
        if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) continue;
        // Get count and sum of valid neighbors
        const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
        float sum= 0.0f;
        if (!iPrecondMode) {
          const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelX ? -iField[x][y][z] : 0.0f);
          const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelY ? -iField[x][y][z] : 0.0f);
          const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (iField[x][y][z]) : (iFieldID == FieldID::IDVelZ ? -iField[x][y][z] : 0.0f);
          if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : iField[x - 1][y][z];
          if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : iField[x + 1][y][z];
          if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : iField[x][y - 1][z];
          if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : iField[x][y + 1][z];
          if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : iField[x][y][z - 1];
          if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : iField[x][y][z + 1];
        }
        // Apply linear expression
        if (iDiffuMode) {
          if (iPrecondMode)
            oField[x][y][z]= 1.0f / (1.0f + diffuVal * (float)count) * iField[x][y][z];            //               [   -D*dt/(h*h)]
          else                                                                                     // [-D*dt/(h*h)] [1+4*D*dt/(h*h)] [-D*dt/(h*h)]
            oField[x][y][z]= (1.0f + diffuVal * (float)count) * iField[x][y][z] - diffuVal * sum;  //               [   -D*dt/(h*h)]
        }
        else {
          if (iPrecondMode)
            oField[x][y][z]= ((voxSize * voxSize) / (float)count) * iField[x][y][z];        //            [-1/(h*h)]
          else                                                                              // [-1/(h*h)] [ 4/(h*h)] [-1/(h*h)]
            oField[x][y][z]= ((float)count * iField[x][y][z] - sum) / (voxSize * voxSize);  //            [-1/(h*h)]
        }
      }
    }
  }
}


// Solve linear system with Conjugate Gradient approach
// References for linear solvers and particularily PCG
// https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// https://services.math.duke.edu/~holee/math361-2020/lectures/Conjugate_gradients.pdf
// https://www3.nd.edu/~zxu2/acms60212-40212-S12/final_project/Linear_solvers_GPU.pdf
// https://github.com/awesson/stable-fluids/tree/master
// https://en.wikipedia.org/wiki/Conjugate_gradient_method
void CompuFluidDyna::ConjugateGradientSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                            const bool iDiffuMode, const float iDiffuCoeff,
                                            const std::vector<std::vector<std::vector<float>>>& iField,
                                            std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[VerboseSolv_____].B()) {
    D.plotLegend.resize(5);
    D.plotLegend[FieldID::IDSmok]= "Diffu S";
    D.plotLegend[FieldID::IDVelX]= "Diffu VX";
    D.plotLegend[FieldID::IDVelY]= "Diffu VY";
    D.plotLegend[FieldID::IDVelZ]= "Diffu VZ";
    D.plotLegend[FieldID::IDPres]= "Proj  P";
    D.plotData.resize(5);
    D.plotData[iFieldID].clear();
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> dField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  dField= rField;
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) printf("\nCG Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("\nCG Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("\nCG Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("\nCG Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("\nCG Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
    D.plotData[iFieldID].push_back(errNew);
  }
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // q = A d
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, dField, qField);
    // alpha = errNew / (d^T q)
    const float denom= ImplicitFieldDotProd(dField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // x = x + alpha d
    ImplicitFieldScale(alpha, dField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ioField= t1Field;
    ApplyBC(iFieldID, ioField);
    // r = r - alpha q
    ImplicitFieldScale(alpha, qField, t0Field);
    ImplicitFieldSub(rField, t0Field, t1Field);
    rField= t1Field;
    // errNew = r^T r
    const float errOld= errNew;
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseSolv_____].B()) {
      printf("%.2e ", errNew);
      D.plotData[iFieldID].push_back(errNew);
    }
    // d = r + (errNew / errOld) * d
    ImplicitFieldScale(errNew / errOld, dField, t0Field);
    ImplicitFieldAdd(rField, t0Field, dField);
  }
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Solve linear system with Gradient descent approach
// Reference on page 55 of https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
void CompuFluidDyna::GradientDescentSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                          const bool iDiffuMode, const float iDiffuCoeff,
                                          const std::vector<std::vector<std::vector<float>>>& iField,
                                          std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[VerboseSolv_____].B()) {
    D.plotLegend.resize(5);
    D.plotLegend[FieldID::IDSmok]= "Diffu S";
    D.plotLegend[FieldID::IDVelX]= "Diffu VX";
    D.plotLegend[FieldID::IDVelY]= "Diffu VY";
    D.plotLegend[FieldID::IDVelZ]= "Diffu VZ";
    D.plotLegend[FieldID::IDPres]= "Proj  P";
    D.plotData.resize(5);
    D.plotData[iFieldID].clear();
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> qField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t1Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) printf("\nGD Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("\nGD Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("\nGD Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("\nGD Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("\nGD Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
    D.plotData[iFieldID].push_back(errNew);
  }
  // Iterate to solve
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= 0.0f) break;
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // q = A r
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, rField, qField);
    // alpha = errNew / (r^T q)
    const float denom= ImplicitFieldDotProd(rField, qField);
    if (denom == 0.0) break;
    const float alpha= errNew / denom;
    // x = x + alpha r
    ImplicitFieldScale(alpha, rField, t0Field);
    ImplicitFieldAdd(ioField, t0Field, t1Field);
    ioField= t1Field;
    ApplyBC(iFieldID, ioField);
    // r = r - alpha q
    ImplicitFieldScale(alpha, qField, t0Field);
    ImplicitFieldSub(rField, t0Field, t1Field);
    rField= t1Field;
    // errNew = r^T r
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseSolv_____].B()) {
      printf("%.2e ", errNew);
      D.plotData[iFieldID].push_back(errNew);
    }
  }
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Solve linear system with an iterative Gauss Seidel scheme
// Solving each equation sequentially by cascading latest solution to next equation
// Run a forward and backward pass in parallel to avoid element ordering bias
// Apply successive overrelaxation coefficient to accelerate convergence
void CompuFluidDyna::GaussSeidelSolve(const int iFieldID, const int iMaxIter, const float iTimeStep,
                                      const bool iDiffuMode, const float iDiffuCoeff,
                                      const std::vector<std::vector<std::vector<float>>>& iField,
                                      std::vector<std::vector<std::vector<float>>>& ioField) {
  // Prepare convergence plot
  if (D.UI[VerboseSolv_____].B()) {
    D.plotLegend.resize(5);
    D.plotLegend[FieldID::IDSmok]= "Diffu S";
    D.plotLegend[FieldID::IDVelX]= "Diffu VX";
    D.plotLegend[FieldID::IDVelY]= "Diffu VY";
    D.plotLegend[FieldID::IDVelZ]= "Diffu VZ";
    D.plotLegend[FieldID::IDPres]= "Proj  P";
    D.plotData.resize(5);
    D.plotData[iFieldID].clear();
  }
  // Allocate fields
  std::vector<std::vector<std::vector<float>>> rField= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<float>>> t0Field= Field::AllocField3D(nX, nY, nZ, 0.0f);
  std::vector<std::vector<std::vector<std::vector<float>>>> FieldT(2);
  // Compute residual error magnitude    r = b - A x    errNew = r · r
  ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
  ApplyBC(iFieldID, t0Field);
  ImplicitFieldSub(iField, t0Field, rField);
  const float errBeg= ImplicitFieldDotProd(rField, rField);
  const float normRHS= ImplicitFieldDotProd(iField, iField);
  float errNew= errBeg;
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) printf("\nGS Diffu S  [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelX) printf("\nGS Diffu VX [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelY) printf("\nGS Diffu VY [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDVelZ) printf("\nGS Diffu VZ [%.2e] ", normRHS);
    if (iFieldID == FieldID::IDPres) printf("\nGS Proj  P  [%.2e] ", normRHS);
    printf("%.2e ", errNew);
    D.plotData[iFieldID].push_back(errNew);
  }
  // Precompute values
  const float diffuVal= iDiffuCoeff * iTimeStep / (voxSize * voxSize);
  const float coeffOverrelax= std::max(D.UI[SolvSOR_________].F(), 0.0f);
  // Iterate to solve with Gauss-Seidel scheme
  for (int idxIter= 0; idxIter < iMaxIter; idxIter++) {
    // Check exit conditions
    if (errNew <= D.UI[SolvTolAbs______].F()) break;
    if (errNew / normRHS <= std::max(D.UI[SolvTolRhs______].F(), 0.0f)) break;
    if (errNew / errBeg <= std::max(D.UI[SolvTolRel______].F(), 0.0f)) break;
    // Initialize fields for forward and backward passes
    FieldT[0]= ioField;
    FieldT[1]= ioField;
    // Execute the two passes in parallel
#pragma omp parallel for
    for (int k= 0; k < 2; k++) {
      // Set the loop settings for the current pass
      const int xBeg= (k == 0) ? 0 : nX - 1;
      const int yBeg= (k == 0) ? 0 : nY - 1;
      const int zBeg= (k == 0) ? 0 : nZ - 1;
      const int xEnd= (k == 0) ? nX : -1;
      const int yEnd= (k == 0) ? nY : -1;
      const int zEnd= (k == 0) ? nZ : -1;
      const int xInc= (k == 0) ? 1 : -1;
      const int yInc= (k == 0) ? 1 : -1;
      const int zInc= (k == 0) ? 1 : -1;
      // Sweep through the field
      for (int x= xBeg; x != xEnd; x+= xInc) {
        for (int y= yBeg; y != yEnd; y+= yInc) {
          for (int z= zBeg; z != zEnd; z+= zInc) {
            // Skip solid or fixed values
            if (Solid[x][y][z]) continue;
            if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
            if (VelBC[x][y][z] && (iFieldID == FieldID::IDVelX || iFieldID == FieldID::IDVelY || iFieldID == FieldID::IDVelZ)) continue;
            if (PreBC[x][y][z] && iFieldID == FieldID::IDPres) continue;
            // Get count and sum of valid neighbors
            const int count= (x > 0) + (y > 0) + (z > 0) + (x < nX - 1) + (y < nY - 1) + (z < nZ - 1);
            float sum= 0.0f;
            const float xBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelX ? -FieldT[k][x][y][z] : 0.0f);
            const float yBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelY ? -FieldT[k][x][y][z] : 0.0f);
            const float zBCVal= (iFieldID == FieldID::IDSmok || iFieldID == FieldID::IDPres) ? (FieldT[k][x][y][z]) : (iFieldID == FieldID::IDVelZ ? -FieldT[k][x][y][z] : 0.0f);
            if (x - 1 >= 0) sum+= Solid[x - 1][y][z] ? xBCVal : FieldT[k][x - 1][y][z];
            if (x + 1 < nX) sum+= Solid[x + 1][y][z] ? xBCVal : FieldT[k][x + 1][y][z];
            if (y - 1 >= 0) sum+= Solid[x][y - 1][z] ? yBCVal : FieldT[k][x][y - 1][z];
            if (y + 1 < nY) sum+= Solid[x][y + 1][z] ? yBCVal : FieldT[k][x][y + 1][z];
            if (z - 1 >= 0) sum+= Solid[x][y][z - 1] ? zBCVal : FieldT[k][x][y][z - 1];
            if (z + 1 < nZ) sum+= Solid[x][y][z + 1] ? zBCVal : FieldT[k][x][y][z + 1];
            // Set new value according to coefficients and flags
            if (count > 0) {
              const float prevVal= FieldT[k][x][y][z];
              if (iDiffuMode) FieldT[k][x][y][z]= (iField[x][y][z] + diffuVal * sum) / (1.0f + diffuVal * (float)count);
              else FieldT[k][x][y][z]= ((voxSize * voxSize) * iField[x][y][z] + sum) / (float)count;
              FieldT[k][x][y][z]= prevVal + coeffOverrelax * (FieldT[k][x][y][z] - prevVal);
            }
          }
        }
      }
    }
    // Recombine forward and backward passes
    for (int x= 0; x < nX; x++)
      for (int y= 0; y < nY; y++)
        for (int z= 0; z < nZ; z++)
          ioField[x][y][z]= (FieldT[0][x][y][z] + FieldT[1][x][y][z]) / 2.0f;
    // Compute residual error magnitude    r = b - A x    errNew = r · r
    ImplicitFieldLaplacianMatMult(iFieldID, iTimeStep, iDiffuMode, iDiffuCoeff, false, ioField, t0Field);
    ApplyBC(iFieldID, t0Field);
    ImplicitFieldSub(iField, t0Field, rField);
    errNew= ImplicitFieldDotProd(rField, rField);
    // Error plot
    if (D.UI[VerboseSolv_____].B()) {
      printf("%.2e ", errNew);
      D.plotData[iFieldID].push_back(errNew);
    }
  }
  // Error plot
  if (D.UI[VerboseSolv_____].B()) {
    if (iFieldID == FieldID::IDSmok) Dum0= rField;
    if (iFieldID == FieldID::IDVelX) Dum1= rField;
    if (iFieldID == FieldID::IDVelY) Dum2= rField;
    if (iFieldID == FieldID::IDVelZ) Dum3= rField;
    if (iFieldID == FieldID::IDPres) Dum4= rField;
  }
}


// Add external forces to velocity field
// vel ⇐ vel + Δt * F / ρ
void CompuFluidDyna::ExternalForces() {
  // Update velocities based on applied external forces
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z] || VelBC[x][y][z]) continue;
        VelZ[x][y][z]+= D.UI[TimeStep________].F() * D.UI[CoeffGravi______].F() * Smok[x][y][z] / fluidDensity;
      }
    }
  }
}


// Project velocity field into a solenoidal/divergence-free field
// 1. Compute RHS based on divergence
// RHS = -(ρ / Δt) × ∇ · vel
// 2. Solve for pressure in pressure Poisson equation
// (-∇²) press = RHS
// 3. Update velocity field by subtracting gradient of pressure
// vel ⇐ vel - (Δt / ρ) × ∇ press
// References for pressure poisson equation and incompressiblity projection
// https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
// https://mycourses.aalto.fi/pluginfile.php/891524/mod_folder/content/0/Lecture03_Pressure.pdf
// https://barbagroup.github.io/essential_skills_RRC/numba/4/#application-pressure-poisson-equation
// http://www.thevisualroom.com/poisson_for_pressure.html
// https://github.com/barbagroup/CFDPython
void CompuFluidDyna::ProjectField(const int iIter, const float iTimeStep,
                                  std::vector<std::vector<std::vector<float>>>& ioVelX,
                                  std::vector<std::vector<std::vector<float>>>& ioVelY,
                                  std::vector<std::vector<std::vector<float>>>& ioVelZ) {
  // Compute divergence for RHS
  ComputeVelocityDivergence();
  // Reset pressure guess to test convergence
  if (D.UI[CoeffProj_______].I() == 2) {
    Pres= Field::AllocField3D(nX, nY, nZ, 0.0f);
    ApplyBC(FieldID::IDPres, Pres);
  }
  // Solve for pressure in the pressure Poisson equation
  if (D.UI[SolvType________].I() == 0) {
    GaussSeidelSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }
  else if (D.UI[SolvType________].I() == 1) {
    GradientDescentSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }
  else {
    ConjugateGradientSolve(FieldID::IDPres, iIter, iTimeStep, false, 0.0f, Dive, Pres);
  }

  // Update velocities based on local pressure gradient
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z] || VelBC[x][y][z]) continue;
        // Subtract pressure gradient to remove divergence
        if (x - 1 >= 0 && !Solid[x - 1][y][z]) ioVelX[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x - 1][y][z]) / (2.0f * voxSize);
        if (y - 1 >= 0 && !Solid[x][y - 1][z]) ioVelY[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y - 1][z]) / (2.0f * voxSize);
        if (z - 1 >= 0 && !Solid[x][y][z - 1]) ioVelZ[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z] - Pres[x][y][z - 1]) / (2.0f * voxSize);
        if (x + 1 < nX && !Solid[x + 1][y][z]) ioVelX[x][y][z]-= iTimeStep / fluidDensity * (Pres[x + 1][y][z] - Pres[x][y][z]) / (2.0f * voxSize);
        if (y + 1 < nY && !Solid[x][y + 1][z]) ioVelY[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y + 1][z] - Pres[x][y][z]) / (2.0f * voxSize);
        if (z + 1 < nZ && !Solid[x][y][z + 1]) ioVelZ[x][y][z]-= iTimeStep / fluidDensity * (Pres[x][y][z + 1] - Pres[x][y][z]) / (2.0f * voxSize);
      }
    }
  }
}


// Trilinearly interpolate the field value at the given position
float CompuFluidDyna::TrilinearInterpolation(const float iPosX, const float iPosY, const float iPosZ,
                                             const std::vector<std::vector<std::vector<float>>>& iFieldRef) {
  // Get floor and ceil voxel indices
  const int x0= std::min(std::max((int)std::floor(iPosX), 0), nX - 1);
  const int y0= std::min(std::max((int)std::floor(iPosY), 0), nY - 1);
  const int z0= std::min(std::max((int)std::floor(iPosZ), 0), nZ - 1);
  const int x1= std::min(std::max((int)std::ceil(iPosX), 0), nX - 1);
  const int y1= std::min(std::max((int)std::ceil(iPosY), 0), nY - 1);
  const int z1= std::min(std::max((int)std::ceil(iPosZ), 0), nZ - 1);
  // Get floor and ceil voxel weights
  const float xWeight1= iPosX - (float)x0;
  const float yWeight1= iPosY - (float)y0;
  const float zWeight1= iPosZ - (float)z0;
  const float xWeight0= 1.0f - xWeight1;
  const float yWeight0= 1.0f - yWeight1;
  const float zWeight0= 1.0f - zWeight1;
  // Compute the weighted sum
  return iFieldRef[x0][y0][z0] * (xWeight0 * yWeight0 * zWeight0) +
         iFieldRef[x0][y0][z1] * (xWeight0 * yWeight0 * zWeight1) +
         iFieldRef[x0][y1][z0] * (xWeight0 * yWeight1 * zWeight0) +
         iFieldRef[x0][y1][z1] * (xWeight0 * yWeight1 * zWeight1) +
         iFieldRef[x1][y0][z0] * (xWeight1 * yWeight0 * zWeight0) +
         iFieldRef[x1][y0][z1] * (xWeight1 * yWeight0 * zWeight1) +
         iFieldRef[x1][y1][z0] * (xWeight1 * yWeight1 * zWeight0) +
         iFieldRef[x1][y1][z1] * (xWeight1 * yWeight1 * zWeight1);
}


// Apply semi-Lagrangian advection along the velocity field
// vel ⇐ vel - Δt (vel · ∇) vel
// smo ⇐ smo - Δt (vel · ∇) smo
// References for MacCormack backtracking scheme
// https://commons.wikimedia.org/wiki/File:Backtracking_maccormack.png
// https://physbam.stanford.edu/~fedkiw/papers/stanford2006-09.pdf
// https://github.com/NiallHornFX/StableFluids3D-GL/blob/master/src/fluidsolver3d.cpp
void CompuFluidDyna::AdvectField(const int iFieldID, const float iTimeStep,
                                 const std::vector<std::vector<std::vector<float>>>& iVelX,
                                 const std::vector<std::vector<std::vector<float>>>& iVelY,
                                 const std::vector<std::vector<std::vector<float>>>& iVelZ,
                                 std::vector<std::vector<std::vector<float>>>& ioField) {
  // Adjust the source field to make solid voxels have a value dependant on their non-solid neighbors
  std::vector<std::vector<std::vector<float>>> sourceField= ioField;
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (!Solid[x][y][z]) continue;
        int count= 0;
        float sum= 0.0f;
        if (x - 1 >= 0 && !Solid[x - 1][y][z] && ++count) sum+= ioField[x - 1][y][z];
        if (y - 1 >= 0 && !Solid[x][y - 1][z] && ++count) sum+= ioField[x][y - 1][z];
        if (z - 1 >= 0 && !Solid[x][y][z - 1] && ++count) sum+= ioField[x][y][z - 1];
        if (x + 1 < nX && !Solid[x + 1][y][z] && ++count) sum+= ioField[x + 1][y][z];
        if (y + 1 < nY && !Solid[x][y + 1][z] && ++count) sum+= ioField[x][y + 1][z];
        if (z + 1 < nZ && !Solid[x][y][z + 1] && ++count) sum+= ioField[x][y][z + 1];
        if (iFieldID == FieldID::IDSmok) sourceField[x][y][z]= (count > 0) ? sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelX) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelY) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
        if (iFieldID == FieldID::IDVelZ) sourceField[x][y][z]= (count > 0) ? -sum / (float)count : 0.0f;
      }
    }
  }
  // Sweep through the field
  for (int x= 0; x < nX; x++) {
#pragma omp parallel for
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        AdvX[x][y][z]= AdvY[x][y][z]= AdvZ[x][y][z]= 0.0f;
        // Skip solid or fixed values
        if (Solid[x][y][z]) continue;
        if (SmoBC[x][y][z] && iFieldID == FieldID::IDSmok) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelX) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelY) continue;
        if (VelBC[x][y][z] && iFieldID == FieldID::IDVelZ) continue;
        // Find source position for active voxel using naive linear backtracking scheme
        const Vec::Vec3<float> posEnd((float)x, (float)y, (float)z);
        const Vec::Vec3<float> velEnd(iVelX[x][y][z], iVelY[x][y][z], iVelZ[x][y][z]);
        Vec::Vec3<float> posBeg= posEnd - iTimeStep * velEnd / voxSize;
        // Iterative source position correction with 2nd order MacCormack scheme
        int correcMaxIter= std::max(D.UI[CoeffAdvec______].I() - 1, 0);
        for (int iter= 0; iter < correcMaxIter; iter++) {
          const float velBegX= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelX);
          const float velBegY= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelY);
          const float velBegZ= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], iVelZ);
          const Vec::Vec3<float> velBeg(velBegX, velBegY, velBegZ);
          const Vec::Vec3<float> vecErr= posEnd - (posBeg + iTimeStep * velBeg / voxSize);
          posBeg= posBeg + vecErr / 2.0f;
        }
        // Save source vector for display
        AdvX[x][y][z]= posBeg[0] - posEnd[0];
        AdvY[x][y][z]= posBeg[1] - posEnd[1];
        AdvZ[x][y][z]= posBeg[2] - posEnd[2];
        // Trilinear interpolation at source position
        ioField[x][y][z]= TrilinearInterpolation(posBeg[0], posBeg[1], posBeg[2], sourceField);
      }
    }
  }
}


// Counteract energy dissipation and introduce turbulent-like behavior by amplifying vorticity on small scales
// https://github.com/awesson/stable-fluids/tree/master
// https://github.com/woeishi/StableFluids/blob/master/StableFluid3d.cpp
// vel ⇐ vel + Δt * TODO write formula
void CompuFluidDyna::VorticityConfinement(const float iTimeStep, const float iVortiCoeff,
                                          std::vector<std::vector<std::vector<float>>>& ioVelX,
                                          std::vector<std::vector<std::vector<float>>>& ioVelY,
                                          std::vector<std::vector<std::vector<float>>>& ioVelZ) {
  // Compute curl and vorticity from the velocity field
  ComputeVelocityCurlVorticity();
  // Amplify non-zero vorticity
  if (iVortiCoeff > 0.0f) {
    for (int x= 0; x < nX; x++) {
      for (int y= 0; y < nY; y++) {
        for (int z= 0; z < nZ; z++) {
          if (Solid[x][y][z] || VelBC[x][y][z]) continue;
          // Gradient of vorticity with zero derivative at solid interface or domain boundary
          Vec::Vec3<float> vortGrad(0.0f, 0.0f, 0.0f);
          if (x - 1 >= 0 && !Solid[x - 1][y][z]) vortGrad[0]+= (Vort[x][y][z] - Vort[x - 1][y][z]) / (2.0f * voxSize);
          if (y - 1 >= 0 && !Solid[x][y - 1][z]) vortGrad[1]+= (Vort[x][y][z] - Vort[x][y - 1][z]) / (2.0f * voxSize);
          if (z - 1 >= 0 && !Solid[x][y][z - 1]) vortGrad[2]+= (Vort[x][y][z] - Vort[x][y][z - 1]) / (2.0f * voxSize);
          if (x + 1 < nX && !Solid[x + 1][y][z]) vortGrad[0]+= (Vort[x + 1][y][z] - Vort[x][y][z]) / (2.0f * voxSize);
          if (y + 1 < nY && !Solid[x][y + 1][z]) vortGrad[1]+= (Vort[x][y + 1][z] - Vort[x][y][z]) / (2.0f * voxSize);
          if (z + 1 < nZ && !Solid[x][y][z + 1]) vortGrad[2]+= (Vort[x][y][z + 1] - Vort[x][y][z]) / (2.0f * voxSize);
          // Amplification of small scale vorticity by following current curl
          if (vortGrad.norm() > 0.0f) {
            const float dVort_dx_scaled= iVortiCoeff * vortGrad[0] / vortGrad.norm();
            const float dVort_dy_scaled= iVortiCoeff * vortGrad[1] / vortGrad.norm();
            const float dVort_dz_scaled= iVortiCoeff * vortGrad[2] / vortGrad.norm();
            ioVelX[x][y][z]+= iTimeStep * (dVort_dy_scaled * CurZ[x][y][z] - dVort_dz_scaled * CurY[x][y][z]);
            ioVelY[x][y][z]+= iTimeStep * (dVort_dz_scaled * CurX[x][y][z] - dVort_dx_scaled * CurZ[x][y][z]);
            ioVelZ[x][y][z]+= iTimeStep * (dVort_dx_scaled * CurY[x][y][z] - dVort_dy_scaled * CurX[x][y][z]);
          }
        }
      }
    }
  }
}


// Compute RHS of pressure poisson equation as negative divergence scaled by density and timestep
// https://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
// RHS = -(ρ / Δt) × ∇ · vel
// TODO implement correction to avoid checkerboard due to odd-even decoupling in pure pressure driven flows
// References for Rhie Chow correction
// https://youtu.be/yqZ59Xn_aF8 Checkerboard oscillations
// https://youtu.be/PmEUiUB8ETk Deriving the correction
// https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2007/rhiechow.pdf OpenFOAM variant
// https://mustafabhotvawala.com/wp-content/uploads/2020/11/MB_rhieChow-1.pdf
void CompuFluidDyna::ComputeVelocityDivergence() {
  // // Precompute pressure gradient for Rhie and Chow correction
  // std::vector<std::vector<std::vector<float>>> PresGradX;
  // std::vector<std::vector<std::vector<float>>> PresGradY;
  // std::vector<std::vector<std::vector<float>>> PresGradZ;
  // if (iUseRhieChow) {
  //   PresGradX= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   PresGradY= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   PresGradZ= Field::AllocField3D(nX, nY, nZ, 0.0f);
  //   for (int x= 0; x < nX; x++) {
  //     for (int y= 0; y < nY; y++) {
  //       for (int z= 0; z < nZ; z++) {
  //         if (Solid[x][y][z]) continue;
  //         // Pressure gradient with zero derivative at solid interface or domain boundary
  //         if (x - 1 >= 0 && !Solid[x - 1][y][z]) PresGradX[x][y][z]+= (Pres[x][y][z] - Pres[x - 1][y][z]) / (2.0f * voxSize);
  //         if (y - 1 >= 0 && !Solid[x][y - 1][z]) PresGradY[x][y][z]+= (Pres[x][y][z] - Pres[x][y - 1][z]) / (2.0f * voxSize);
  //         if (z - 1 >= 0 && !Solid[x][y][z - 1]) PresGradZ[x][y][z]+= (Pres[x][y][z] - Pres[x][y][z - 1]) / (2.0f * voxSize);
  //         if (x + 1 < nX && !Solid[x + 1][y][z]) PresGradX[x][y][z]+= (Pres[x + 1][y][z] - Pres[x][y][z]) / (2.0f * voxSize);
  //         if (y + 1 < nY && !Solid[x][y + 1][z]) PresGradY[x][y][z]+= (Pres[x][y + 1][z] - Pres[x][y][z]) / (2.0f * voxSize);
  //         if (z + 1 < nZ && !Solid[x][y][z + 1]) PresGradZ[x][y][z]+= (Pres[x][y][z + 1] - Pres[x][y][z]) / (2.0f * voxSize);
  //       }
  //     }
  //   }
  // }
  // Compute divergence of velocity field
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        if (Solid[x][y][z]) Dive[x][y][z]= 0.0f;
        if (PreBC[x][y][z]) Dive[x][y][z]= PresForced[x][y][z];
        if (Solid[x][y][z] || PreBC[x][y][z]) continue;
        // Classical linear interpolation for face velocities with same velocity at domain boundary and zero velocity at solid interface
        float velXN= (x - 1 >= 0) ? ((Solid[x - 1][y][z]) ? (0.0f) : ((VelX[x][y][z] + VelX[x - 1][y][z]) / 2.0f)) : (VelX[x][y][z]);
        float velYN= (y - 1 >= 0) ? ((Solid[x][y - 1][z]) ? (0.0f) : ((VelY[x][y][z] + VelY[x][y - 1][z]) / 2.0f)) : (VelY[x][y][z]);
        float velZN= (z - 1 >= 0) ? ((Solid[x][y][z - 1]) ? (0.0f) : ((VelZ[x][y][z] + VelZ[x][y][z - 1]) / 2.0f)) : (VelZ[x][y][z]);
        float velXP= (x + 1 < nX) ? ((Solid[x + 1][y][z]) ? (0.0f) : ((VelX[x + 1][y][z] + VelX[x][y][z]) / 2.0f)) : (VelX[x][y][z]);
        float velYP= (y + 1 < nY) ? ((Solid[x][y + 1][z]) ? (0.0f) : ((VelY[x][y + 1][z] + VelY[x][y][z]) / 2.0f)) : (VelY[x][y][z]);
        float velZP= (z + 1 < nZ) ? ((Solid[x][y][z + 1]) ? (0.0f) : ((VelZ[x][y][z + 1] + VelZ[x][y][z]) / 2.0f)) : (VelZ[x][y][z]);
        // // Rhie and Chow correction terms
        // if (iUseRhieChow) {
        //   // Subtract pressure gradients with neighboring cells
        //   velXN-= D.UI[CoeffProj1______].F() * ((x - 1 >= 0 && !Solid[x - 1][y][z]) ? ((Pres[x][y][z] - Pres[x - 1][y][z]) / voxSize) : (0.0f));
        //   velYN-= D.UI[CoeffProj1______].F() * ((y - 1 >= 0 && !Solid[x][y - 1][z]) ? ((Pres[x][y][z] - Pres[x][y - 1][z]) / voxSize) : (0.0f));
        //   velZN-= D.UI[CoeffProj1______].F() * ((z - 1 >= 0 && !Solid[x][y][z - 1]) ? ((Pres[x][y][z] - Pres[x][y][z - 1]) / voxSize) : (0.0f));
        //   velXP-= D.UI[CoeffProj1______].F() * ((x + 1 < nX && !Solid[x + 1][y][z]) ? ((Pres[x + 1][y][z] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   velYP-= D.UI[CoeffProj1______].F() * ((y + 1 < nY && !Solid[x][y + 1][z]) ? ((Pres[x][y + 1][z] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   velZP-= D.UI[CoeffProj1______].F() * ((z + 1 < nZ && !Solid[x][y][z + 1]) ? ((Pres[x][y][z + 1] - Pres[x][y][z]) / voxSize) : (0.0f));
        //   // Add Linear interpolations of pressure gradients with neighboring cells
        //   velXN+= D.UI[CoeffProj2______].F() * ((x - 1 >= 0) ? ((PresGradX[x][y][z] + PresGradX[x - 1][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velYN+= D.UI[CoeffProj2______].F() * ((y - 1 >= 0) ? ((PresGradY[x][y][z] + PresGradY[x][y - 1][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velZN+= D.UI[CoeffProj2______].F() * ((z - 1 >= 0) ? ((PresGradZ[x][y][z] + PresGradZ[x][y][z - 1]) / 2.0f) : (PresGradX[x][y][z]));
        //   velXP+= D.UI[CoeffProj2______].F() * ((x + 1 < nX) ? ((PresGradX[x + 1][y][z] + PresGradX[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velYP+= D.UI[CoeffProj2______].F() * ((y + 1 < nY) ? ((PresGradY[x][y + 1][z] + PresGradY[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        //   velZP+= D.UI[CoeffProj2______].F() * ((z + 1 < nZ) ? ((PresGradZ[x][y][z + 1] + PresGradZ[x][y][z]) / 2.0f) : (PresGradX[x][y][z]));
        // }
        // Divergence based on face velocities scaled by density and timestep  (negated RHS and linear system to have positive diag coeffs)
        Dive[x][y][z]= -fluidDensity / D.UI[TimeStep________].F() * ((velXP - velXN) + (velYP - velYN) + (velZP - velZN)) / voxSize;
      }
    }
  }
}


// Compute curl and vorticity of current velocity field
// curl= ∇ ⨯ vel
// vort= ‖curl‖₂
void CompuFluidDyna::ComputeVelocityCurlVorticity() {
  for (int x= 0; x < nX; x++) {
    for (int y= 0; y < nY; y++) {
      for (int z= 0; z < nZ; z++) {
        CurX[x][y][z]= CurY[x][y][z]= CurZ[x][y][z]= Vort[x][y][z]= 0.0f;
        if (Solid[x][y][z]) continue;
        // Compute velocity cross derivatives considering BC at interface with solid
        float dVely_dx= 0.0f, dVelz_dx= 0.0f, dVelx_dy= 0.0f, dVelz_dy= 0.0f, dVelx_dz= 0.0f, dVely_dz= 0.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVely_dx= ((Solid[x + 1][y][z] ? VelY[x][y][z] : VelY[x + 1][y][z]) - (Solid[x - 1][y][z] ? VelY[x][y][z] : VelY[x - 1][y][z])) / 2.0f;
        if (x - 1 >= 0 && x + 1 < nX) dVelz_dx= ((Solid[x + 1][y][z] ? VelZ[x][y][z] : VelZ[x + 1][y][z]) - (Solid[x - 1][y][z] ? VelZ[x][y][z] : VelZ[x - 1][y][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVelx_dy= ((Solid[x][y + 1][z] ? VelX[x][y][z] : VelX[x][y + 1][z]) - (Solid[x][y - 1][z] ? VelX[x][y][z] : VelX[x][y - 1][z])) / 2.0f;
        if (y - 1 >= 0 && y + 1 < nY) dVelz_dy= ((Solid[x][y + 1][z] ? VelZ[x][y][z] : VelZ[x][y + 1][z]) - (Solid[x][y - 1][z] ? VelZ[x][y][z] : VelZ[x][y - 1][z])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVelx_dz= ((Solid[x][y][z + 1] ? VelX[x][y][z] : VelX[x][y][z + 1]) - (Solid[x][y][z - 1] ? VelX[x][y][z] : VelX[x][y][z - 1])) / 2.0f;
        if (z - 1 >= 0 && z + 1 < nZ) dVely_dz= ((Solid[x][y][z + 1] ? VelY[x][y][z] : VelY[x][y][z + 1]) - (Solid[x][y][z - 1] ? VelY[x][y][z] : VelY[x][y][z - 1])) / 2.0f;
        // Deduce curl and vorticity
        CurX[x][y][z]= dVelz_dy - dVely_dz;
        CurY[x][y][z]= dVelx_dz - dVelz_dx;
        CurZ[x][y][z]= dVely_dx - dVelx_dy;
        Vort[x][y][z]= std::sqrt(CurX[x][y][z] * CurX[x][y][z] + CurY[x][y][z] * CurY[x][y][z] + CurZ[x][y][z] * CurZ[x][y][z]);
      }
    }
  }
}
