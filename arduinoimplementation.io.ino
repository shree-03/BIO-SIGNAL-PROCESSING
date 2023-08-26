#include <SPI.h>
#include<SD.h>
#define PI 3.1415926535897932384626433832795
/*https://archive.physionet.org/cgi-bin/atm/ATM */

void dft(const double inrealn[],double outreal[],double outimg[],double magi[],size_t n);//for n number of samples
float maxi(float a, float b)

double inreal[1500]={0.2593,0.4559,0.394,0.4186,0.4035,0.4127,0.4098,0.496,0.5763,0.6014,0.5985,0.5585,0.5038,0.4934,0.4846,0.467,0.4443,0.4228,0.4139,0.4125,0.4185,0.4184,0.4154,0.4679,0.5582,0.5919,0.5824,0.5419,0.4871,0.4669,0.4529,0.4303,0.417,0.4118,0.4143,0.4171,0.4196,0.4172,0.4109,0.4547,0.5378,0.5769,0.5749,0.533,0.4845,0.4633,0.4457,0.4253,0.4153,0.4088,0.4083,0.4082,0.4125,0.407,0.4155,0.4928,0.5788,0.6097,0.6053,0.5637,0.5123,0.4948,0.4765,0.4599,0.4469,0.4229,0.4092,0.4112,0.4105,0.4112,0.4118,0.4754,0.5696,0.5986,0.5909,0.5435,0.4837,0.4665,0.4595,0.4401,0.4176,0.4083,0.4035,0.3963,0.3884,0.3849,0.3794,0.4233,0.5141,0.5694,0.5773,0.5399,0.4892,0.4731,0.4625,0.4493,0.4262,0.4089,0.4029,0.3992,0.3976,0.3863,0.3804,0.4313,0.5265,0.5699,0.5741,0.5394,0.4932,0.4791,0.4667,0.4417,0.4233,0.4133,0.4119,0.41,0.4124,0.4229,0.4247,0.4955,0.5874,0.6195,0.618,0.5791,0.5215,0.4981,0.4883,0.4675,0.4482,0.4277,0.4137,0.4073,0.4136,0.4146,0.4112,0.4549,0.5446,0.5871,0.587,0.5486,0.4972,0.4805,0.4663,0.4444,0.424,0.413,0.4148,0.4221,0.4231,0.4156,0.4095,0.4526,0.5322,0.5758,0.5771,0.5308,0.4823,0.4664,0.4546,0.4311,0.4191,0.4133,0.4178,0.4227,0.4242,0.4157,0.433,0.5245,0.5924,0.6108,0.6027,0.5569,0.5009,0.4783,0.4629,0.4498,0.4279,0.4119,0.4057,0.4065,0.4062,0.3997,0.4163,0.4995,0.5754,0.5853,0.5753,0.5234,0.4714,0.4541,0.4381,0.4209,0.4093,0.3964,0.3904,0.3926,0.387,0.3789,0.3887,0.4691,0.5481,0.5715,0.5567,0.5124,0.4645,0.4572,0.45,0.4275,0.4094,0.408,0.4069,0.4053,0.3982,0.3889,0.4028,0.471,0.5509,0.5646,0.5334,0.4855,0.4695,0.4591,0.4434,0.4316,0.4208,0.4134,0.4111,0.412,0.4181,0.4177,0.4512,0.5516,0.6028,0.6047,0.5844,0.5364,0.4989,0.4815,0.4628,0.4442,0.4256,0.4156,0.4166,0.4175,0.4162,0.4114,0.4454,0.526,0.5719,0.5789,0.5568,0.5043,0.4768,0.4636,0.4384,0.4198,0.4086,0.3998,0.4036,0.4024,0.3992,0.3924,0.4206,0.5039,0.5614,0.5726,0.543,0.492,0.4746,0.4608,0.4394,0.419,0.4082,0.4013,0.4011,0.3975,0.3954,0.395,0.4598,0.5579,0.5858,0.5799,0.5516,0.5031,0.4853,0.4744,0.4561,0.4337,0.4202,0.4183,0.4187,0.4171,0.411,0.4246,0.5062,0.5816,0.6029,0.5985,0.5546,0.5043,0.4875,0.4714,0.4558,0.4349,0.4178,0.4193,0.4116,0.4034,0.4022,0.4118,0.4855,0.5655,0.5778,0.5706,0.5191,0.4724,0.4593,0.4419,0.4229,0.4093,0.4041,0.4011,0.4014,0.4063,0.398,0.4054,0.4869,0.5836,0.6068,0.5839,0.5339,0.5021,0.5002,0.4998,0.5001,0.4875,0.4656,0.4534,0.4432,0.4236,0.4123,0.4153,0.4736,0.5387,0.5511,0.5192,0.4556,0.4194,0.4106,0.4074,0.401,0.3887,0.3833,0.3863,0.3886,0.389,0.3848,0.4255,0.5331,0.5985,0.6027,0.587,0.5408,0.4941,0.4803,0.4605,0.4414,0.4227,0.416,0.4165,0.414,0.413,0.41,0.4538,0.5466,0.5988,0.6094,0.5871,0.5403,0.5059,0.4963,0.473,0.4562,0.4363,0.423,0.424,0.4216,0.4227,0.4209,0.4452,0.5283,0.5861,0.5921,0.5714,0.5178,0.4814,0.4666,0.4529,0.4322,0.4215,0.4144,0.4118,0.4078,0.4008,0.39,0.4046,0.4776,0.5448,0.5637,0.5449,0.49,0.4617,0.452,0.4313,0.4193,0.4086,0.3997,0.3954,0.3879,0.3831,0.3691,0.3941,0.4744,0.5453,0.565,0.5362,0.496,0.4755,0.4653,0.4541,0.4464,0.4371,0.4289,0.4284,0.4282,0.4273,0.4247,0.4634,0.5516,0.593,0.595,0.571,0.5257,0.4999,0.4828,0.4649,0.4489,0.4402,0.4329,0.431,0.4294,0.424,0.4136,0.4434,0.5299,0.5764,0.5807,0.5631,0.5065,0.478,0.4698,0.4494,0.4394,0.4276,0.4192,0.4146,0.4137,0.4123,0.4056,0.4291,0.5101,0.5705,0.5768,0.5381,0.4933,0.4681,0.4479,0.4288,0.4155,0.4061,0.393,0.389,0.3874,0.3879,0.3964,0.4709,0.586,0.6214,0.6099,0.5722,0.516,0.4937,0.4747,0.4517,0.4333,0.4199,0.4119,0.4097,0.4084,0.4036,0.4095,0.4808,0.5899,0.6348,0.6176,0.5767,0.5118,0.4844,0.4594,0.4358,0.4222,0.4158,0.4123,0.4096,0.4069,0.3995,0.3964,0.4446,0.5513,0.6004,0.5988,0.5475,0.4857,0.4621,0.4482,0.4336,0.4158,0.4017,0.3928,0.3888,0.3848,0.3767,0.3695,0.4269,0.5233,0.56,0.5537,0.5105,0.4668,0.4571,0.4368,0.4192,0.4095,0.405,0.4025,0.4083,0.4068,0.4034,0.4069,0.4887,0.598,0.6314,0.6175,0.574,0.5101,0.4909,0.477,0.453,0.4274,0.4173,0.4115,0.4063,0.4009,0.3867,0.3829,0.4452,0.5531,0.5972,0.5904,0.5487,0.4909,0.4672,0.4529,0.4295,0.4136,0.4015,0.3928,0.3953,0.3949,0.3912,0.3875,0.4318,0.5281,0.5781,0.5783,0.5308,0.4729,0.4586,0.4441,0.4223,0.41,0.3991,0.3911,0.3831,0.3859,0.3857,0.3894,0.4801,0.5772,0.5976,0.5879,0.5427,0.4974,0.4821,0.4666,0.4473,0.4302,0.4179,0.4124,0.4089,0.4127,0.4161,0.4226,0.5077,0.5908,0.6095,0.5944,0.5466,0.4996,0.4864,0.4701,0.4475,0.4275,0.4206,0.416,0.405,0.4034,0.4033,0.4057,0.4737,0.5729,0.5965,0.5882,0.5381,0.4862,0.4709,0.4471,0.4267,0.412,0.4062,0.3991,0.3983,0.4029,0.3959,0.3879,0.4486,0.5476,0.5827,0.5699,0.5204,0.4708,0.4555,0.4384,0.4177,0.4111,0.401,0.3919,0.3946,0.4009,0.3939,0.4127,0.5079,0.6003,0.6215,0.6019,0.5468,0.5037,0.4953,0.4743,0.4507,0.429,0.4201,0.4165,0.4172,0.4188,0.4128,0.4207,0.5043,0.6016,0.6175,0.5988,0.5506,0.4941,0.4758,0.4582,0.4306,0.4109,0.4011,0.3985,0.401,0.4001,0.3955,0.3912,0.4549,0.5546,0.593,0.5863,0.5347,0.4842,0.4673,0.4456,0.42,0.4057,0.4036,0.3952,0.3974,0.3952,0.3882,0.3828,0.4275,0.5305,0.5851,0.5746,0.5167,0.4662,0.4495,0.4345,0.4199,0.4097,0.4031,0.4033,0.4077,0.4069,0.405,0.4072,0.4776,0.589,0.6317,0.6274,0.5792,0.5189,0.4836,0.4754,0.4554,0.4314,0.4115,0.4038,0.4026,0.4046,0.4065,0.3985,0.4291,0.533,0.5926,0.5939,0.5665,0.5089,0.4781,0.4649,0.4458,0.4176,0.4022,0.398,0.3996,0.3976,0.3978,0.3959,0.4067,0.4968,0.5761,0.5832,0.5575,0.4963,0.4594,0.4516,0.4278,0.4081,0.394,0.393,0.3934,0.3919,0.3923,0.3825,0.3877,0.4641,0.5568,0.5765,0.5531,0.4904,0.4523,0.4391,0.415,0.4023,0.3939,0.3943,0.3933,0.397,0.3982,0.3888,0.4075,0.5218,0.6247,0.6413,0.623,0.5656,0.507,0.4933,0.4737,0.4543,0.428,0.4099,0.4067,0.4067,0.4051,0.4019,0.4029,0.4684,0.5834,0.6216,0.6108,0.5658,0.5043,0.4901,0.4798,0.4529,0.4199,0.4064,0.3965,0.388,0.3953,0.3955,0.3818,0.418,0.5279,0.598,0.6063,0.5706,0.4991,0.4679,0.4523,0.4265,0.4067,0.3888,0.3805,0.3805,0.389,0.387,0.3741,0.3808,0.4597,0.5586,0.5824,0.5518,0.492,0.4568,0.4398,0.4167,0.3949,0.3803,0.373,0.3752,0.3814,0.38,0.3761,0.3894,0.4931,0.5888,0.6048,0.5914,0.5393,0.4961,0.4829,0.4616,0.4401,0.4205,0.4025,0.3989,0.3904,0.3934,0.3987,0.3945,0.4458,0.5631,0.6046,0.5964,0.5574,0.5005,0.4822,0.4679,0.4443,0.4233,0.4042,0.3913,0.3869,0.3895,0.4026,0.4034,0.4075,0.5021,0.5926,0.6058,0.5843,0.5215,0.4771,0.4707,0.4493,0.4245,0.4059,0.3948,0.3943,0.3958,0.399,0.3928,0.3833,0.4411,0.5477,0.5928,0.584,0.5306,0.4774,0.4586,0.441,0.4192,0.4001,0.388,0.379,0.3849,0.3904,0.3943,0.3938,0.4471,0.5567,0.6146,0.6153,0.5843,0.524,0.499,0.4836,0.4669,0.4504,0.4221,0.411,0.4072,0.405,0.4051,0.4005,0.4416,0.5663,0.6244,0.6189,0.5798,0.5206,0.4988,0.4876,0.4676,0.4429,0.4235,0.4063,0.4044,0.3998,0.3999,0.4024,0.4027,0.4724,0.5769,0.6064,0.5931,0.5396,0.4919,0.4841,0.4663,0.4458,0.4167,0.4046,0.3981,0.4007,0.4054,0.3959,0.383,0.4065,0.5076,0.5839,0.589,0.5587,0.4976,0.473,0.4665,0.4412,0.4234,0.4034,0.3983,0.3982,0.3926,0.3964,0.3834,0.3912,0.4798,0.56,0.5739,0.5542,0.5032,0.483,0.4788,0.4587,0.4405,0.4263,0.4126,0.4133,0.4122,0.413,0.4076,0.4044,0.4826,0.5807,0.6128,0.6066,0.5638,0.5129,0.4981,0.4935,0.4747,0.452,0.43,0.4155,0.4134,0.4118,0.4131,0.4058,0.4337,0.5391,0.5936,0.5941,0.5787,0.5312,0.4965,0.4782,0.4626,0.4479,0.4263,0.4107,0.4122,0.412,0.4092,0.4072,0.4139,0.4857,0.5709,0.593,0.5688,0.5079,0.4742,0.4675,0.4537,0.4386,0.4206,0.4054,0.3985,0.4039,0.405,0.3981,0.3893,0.4455,0.5487,0.5876,0.5733,0.5187,0.4728,0.4535,0.4392,0.4186,0.4064,0.399,0.3965,0.4022,0.404,0.4059,0.4029,0.4428,0.5473,0.6019,0.5992,0.5753,0.5202,0.4896,0.4867,0.4721,0.453,0.4305,0.4122,0.4058,0.4046,0.4028,0.4028,0.4182,0.5164,0.5886,0.5882,0.5723,0.5202,0.4918,0.483,0.4668,0.4541,0.4337,0.4151,0.4104,0.4008,0.3883,0.3886,0.3919,0.3843,0.3852,0.3825,0.3819,0.3865,0.3842,0.3895,0.3977,0.4003,0.4044,0.4087,0.4065,0.4072,0.4072,0.4092,0.4115,0.5088,0.6399,0.6749,0.668,0.6188,0.5667,0.5392,0.5125,0.4883,0.4674,0.4462,0.4252,0.4087,0.4072,0.4011,0.3999,0.4063,0.4656,0.5642,0.5969,0.5901,0.5528,0.501,0.4889,0.4769,0.4562,0.4433,0.4349,0.4265,0.4168,0.4108,0.4105,0.4093,0.4398,0.5116,0.5533,0.5491,0.532,0.4936,0.4744,0.4714,0.4601,0.4458,0.4316,0.4224,0.4134,0.4097,0.4111,0.4054,0.4114,0.4671,0.528,0.5441,0.5354,0.4965,0.4712,0.4638,0.4557,0.4379,0.4242,0.4159,0.4104,0.4123,0.4102,0.4103,0.4112,0.4342,0.504,0.5439,0.5392,0.5154,0.4804,0.4686,0.4584,0.4417,0.4263,0.4104,0.399,0.3939,0.3966,0.3885,0.377,0.402,0.5044,0.5822,0.5918,0.5711,0.5192,0.4874,0.4792,0.4638,0.4477,0.4354,0.4151,0.4034,0.4015,0.4024,0.3991,0.4106,0.4913,0.5891,0.6105,0.5974,0.555,0.5028,0.4885,0.4802,0.4623,0.4376,0.4243,0.4103,0.4025,0.4027,0.3954,0.3884,0.4355,0.5417,0.5899,0.5915,0.5609,0.4994,0.4777,0.469,0.4517,0.4282,0.4094,0.4009,0.4103,0.4139,0.4138,0.404,0.4062,0.4826,0.5705,0.5902,0.5789,0.5282,0.493,0.4817,0.4672,0.4525,0.4376,0.419,0.4144,0.4123,0.4131,0.411,0.4066,0.4524,0.5424,0.5861,0.5856,0.5406,0.4868,0.4715,0.4596,0.45,0.4279,0.4175,0.415,0.4132,0.4147,0.414,0.4089,0.4375,0.5384,0.6016,0.6002,0.5817,0.5366,0.4999,0.4945,0.4732,0.4582,0.443,0.4198,0.4154,0.4105,0.4179,0.4209,0.4185,0.4752,0.5572,0.5874,0.5769,0.5357,0.4961,0.4885,0.4721,0.4588,0.4419,0.4243,0.4168,0.4179,0.4162,0.4149,0.4104,0.4306,0.5136,0.5705,0.5767,0.5594,0.516,0.485,0.4774,0.4616,0.4475,0.43,0.4229,0.4216,0.4215,0.4231,0.4212,0.4192,0.4618,0.5484,0.5784,0.5695,0.5332,0.4884,0.4738,0.4616,0.4406,0.4231,0.4229,0.4227,0.423,0.4218,0.4227,0.4175,0.4472,0.5336,0.5907,0.5969,0.5733,0.5212,0.4942,0.4834,0.4717,0.4631,0.4422,0.4259,0.4207,0.4184,0.4204,0.4198,0.4183,0.4665,0.5593,0.5877,0.5821,0.5462,0.5022,0.4887,0.4766,0.4689,0.4495,0.4264,0.4201,0.4199,0.4176,0.4161,0.4158,0.4283,0.4994,0.5712,0.5834,0.5657,0.5174,0.4891,0.4785,0.4621,0.4488,0.429,0.416,0.4172,0.4122,0.4171,0.4181,0.4096,0.4372,0.5148,0.5684,0.5762,0.5581,0.5042,0.4794,0.4795,0.4674,0.4488,0.4341,0.4198,0.4172,0.4152,0.4119,0.4164,0.4105,0.456,0.5367,0.5679,0.5629,0.5251,0.4849,0.4761,0.463,0.4443,0.4302,0.4167,0.4146,0.4181,0.416,0.4166,0.4101,0.425,0.4984,0.5739,0.5875,0.5771,0.5365,0.4966,0.4861,0.477,0.4638,0.4473,0.4299,0.4182,0.4187,0.418,0.4164,0.4145,0.4305,0.5019,0.5706,0.5825,0.5629,0.5167,0.4863,0.4764,0.4668,0.4476,0.4336,0.4184,0.4156,0.4178,0.418,0.4155,0.4108,0.4484,0.5273,0.5693,0.5746,0.5464,0.4993,0.4813,0.4785,0.4585,0.4374,0.4246,0.4156,0.4179,0.4194,0.4157,0.4144,0.4145,0.4691,0.5498,0.5759,0.5695,0.5313,0.4939,0.4831,0.4773,0.4624,0.4426,0.4232,0.4177,0.4179,0.4208,0.4129,0.4105,0.4115,0.5142,0.5339,0.6241
};

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
}

void loop() {
  // preprocessing 
  //mean 
  int n=1500;

  /*for(int i=0;i<n;i++){
    Serial.print(i);
    Serial.print(" ");
    Serial.println(inreal[i]);
  }// for plotting the original signal*/
  
  double sum1=0;
  for(int j=0;j<n;j++){ 
    sum1=sum1+inreal[j];
    
  }
  //Serial.println(sum1);
  
  //find mean
  double mean=sum1/n;
  //Serial.println(mean);          //all clear
  
  for(int i=0;i<n;i++)
  {
    inreal[i]=inreal[i]-mean;
  }
  /*for(int i=0;i<n;i++)
  {
     Serial.println( inreal[i]);
  }*/
  
  //find max value
  double maxVal = inreal[0];
  
  for (int i = 1; i < n; i++) {
      if (inreal[i] > maxVal) {
         maxVal = inreal[i];
      }
      else{
        //
      }
   }
   //Serial.println(maxVal);  //fine
   
   //normalise the signal
   for(int j=0;j<n;j++){
     inreal[j]=inreal[j]/maxVal;
     //Serial.println(inreal[j]);  //fine
   }

   //this is for coefficient calculation for resonator
   double a1=-2*0.997*cos(0.3*2*PI);
   double a2=0.997*0.997;
   double b=(1-0.997)*sqrt(1+(0.997*0.997)-2*0.997*cos(2*0.3));
   

  //resonator
  double inrealn[1500]={};
  
  for(int j=2;j<n;j++){
    inrealn[j]=-a1*inrealn[j-1]-a2*inrealn[j-2]+b*inreal[j];
  }
  
  //dft
  //double inmg[250]={};
  double outreal[1500]={};
  double outimg[1500]={};
  double magi[1500]={};
  
  dft(inrealn,outreal,outimg,magi,n);//works good
  //Serial.println(magi[],6);  
  /*for(int i=0;i<n;i++){
    //Serial.print(i);
    //Serial.print(",");
    Serial.println(magi[i],6);   //
  }// this is for plotting magnitude spectrum*/

  //find maximum peak in spectrum
  double fmaxi = -100000.00;
  
  for (int i = 0; i < n; i++) {
    fmaxi = maxi(magi[i],fmaxi);
    //Serial.println(fmaxi);
    
    // delay(10);
  }
  
  
  //find breathing rate
  double resp_rate=fmaxi*60;
  Serial.println(resp_rate);
}

//for dft
void dft(const double inrealn[],double outreal[],double outimg[],double magi[],size_t n){//double inmg is removed ass it is asumed to be 0
  
    for (int k=0;k<n;k++){
      double sumreal=0;
      double sumimg=0;
      for (int t=0;t<n;t++){
          double angles=2*PI*t*k/n;
          sumreal+=inrealn[t]*cos(angles);
          sumimg+=-inrealn[t]*sin(angles);
         
      }
      
      outreal[k]=sumreal;
      outimg[k]=sumimg;
      magi[k]=sqrt(pow(outreal[k],2)+pow(outimg[k],2));
      //Serial.println(k); 
      //Serial.println(','); 
      //Serial.println(magi[k],6);  //works good
    }
}

float maxi(float a, float b){
    if(a>b){
        return a;
    }
    else{
        return b;
    }
    
}