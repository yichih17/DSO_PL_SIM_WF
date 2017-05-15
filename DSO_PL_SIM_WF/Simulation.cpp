#include"define.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<random>
#include<time.h>

#define outputPAT 1									//1: output PAT to txt file ; 0: store in program memory
#define outputUEinfo 0								//1: output UE information to txt file ; 0: store in program memory

using namespace std;

int WIFICQIRange[] = { 185, 152, 133, 109, 84, 64, 60, 56 };
double ap_capacity[8] = { 6500, 13000, 19500, 26000, 39000, 52000, 58500, 65000 };	//65Mbps = 65000000bps = 65000 bits/ms
UE UEList[UEnumber];
int DB50_UEnumber = 0;
int DB100_UEnumber = 0;
int DB300_UEnumber = 0;
void OutputResult(string Scheme, SimulationResult *Result);


template <class T>
void uniformdistribution(T* equip)
{
	std::random_device rd;							//integer random number generator that produces non-deterministic random numbers. 
	std::mt19937 gen(rd());							//Mersenne Twister 19937 generator, generate a random number seed
	std::uniform_real_distribution<> theta(0, 360);	//definition of a uniform distribution range, a random number between 0 and 360
	std::uniform_real_distribution<> k(0, 1);		//definition of a uniform distribution range, a random number between 0 and 1
	double r = radius_AP * sqrt(k(gen));			//random a angle and random a radius, to gennerate a coordinate for UE
	double angel = (theta(gen));
	equip->coor_X = r * std::sin(angel);
	equip->coor_Y = r * std::cos(angel);
}

double exponentially_Distributed(double x)
{
	double y, z;
	y = (double)(rand() + 1) / (double)(RAND_MAX + 1);
	z = (double)log(y) * (double)(-1 / x);
	return z;
}

double getDistance(double coor_x, double coor_y) { return sqrt(pow(coor_x, 2) + pow(coor_y, 2)); }

int getCQI(UE *u)
{
	double distance = getDistance(u->coor_X, u->coor_Y);
	int CQI = 0;
	for (int i = 0; i < 8; i++)
	{
		if (distance <= WIFICQIRange[i])
			CQI++;
		else
			break;
	}
	return CQI;
}

// int轉string
string IntToString(int &i)
{
	string s;
	stringstream ss(s);
	ss << i;
	return ss.str();
}

void Simulation_Result(UE *UEList, SimulationResult *Result)
{
	//Queueing model calculation
	double AvgSystemTime = 0.0;
	double Xj = 0.0;
	double Xj2 = 0.0;
	double Xj_paper = 0.0;
	double Xj2_paper = 0.0;
	double lambda = 0.0;
	for (int i = 0; i < UEnumber; i++)
		lambda = lambda + UEList[i].lambdai;
	for (int i = 0; i < UEnumber; i++)
	{
		double weight_i = UEList[i].lambdai / lambda;

		double Xij = UEList[i].packet_size / (ap_capacity[UEList[i].CQI - 1]);
		Xj = Xj + Xij * weight_i;
		Xj2 = Xj2 + pow(Xij, 2) * weight_i;

		double Xij_paper = UEList[i].packet_size / (ap_capacity[UEList[i].CQI - 1] / UEnumber);
		Xj_paper = Xj_paper + Xij_paper * weight_i;
		Xj2_paper = Xj2_paper + pow(Xij_paper, 2) * weight_i;
	}
	//double rho = lambda * Xj;
	Result->AvgSystemTime_paper = Xj_paper + lambda * Xj2_paper / (1 - lambda * Xj_paper);
	double rho = lambda * Xj;
	double right = lambda * Xj2 / (1 - lambda * Xj);
	Result->AvgSystemTime = Xj + lambda * Xj2 / (1 - lambda * Xj);

	// 計算整體的throughput、delay、schedule packet數、discard packet數
	double DelayTemp = 0.0;
	double SystemTimeTemp = 0.0;
	double Type1_DelayTemp = 0.0;
	double Type2_DelayTemp = 0.0;
	double Type3_DelayTemp = 0.0;
	for (int i = 0; i<UEnumber; i++)
	{
		Result->TotalThroughput = Result->TotalThroughput + Result->Throughput[i];
		DelayTemp += Result->Delay[i];
		SystemTimeTemp += Result->SystemTime[i];
		Result->TotalSchedulePacketNum = Result->TotalSchedulePacketNum + Result->SchedulePackerNum[i];
		Result->TotalDiscardPacketNum = Result->TotalDiscardPacketNum + Result->DiscardPacketNum[i];
	}
	Result->AverageThroughput = Result->TotalThroughput / UEnumber;
	Result->AverageDelay = DelayTemp / Result->TotalSchedulePacketNum;
	Result->AverageSystemTime = SystemTimeTemp / Result->TotalSchedulePacketNum;
	Result->PacketLossRatio = ((double)Result->TotalDiscardPacketNum / (double)(Result->TotalSchedulePacketNum + Result->TotalDiscardPacketNum)) * 100;

	//// 計算typ1(VoIP)的throughput、delay、schedule packet數、discard packet數、rate滿意度、delay滿意度
	//if (DB50_UEnumber > 0)
	//{
	//	for (int i = 0; i<DB50_UEnumber; i++)
	//	{
	//		Result->Type1_TotalThroughput = Result->Type1_TotalThroughput + Result->Throughput[i];
	//		Type1_DelayTemp = Type1_DelayTemp + (Result->Delay[i] / Result->SchedulePackerNum[i]);
	//		Result->Type1_SchedulePacketNum = Result->Type1_SchedulePacketNum + Result->SchedulePackerNum[i];
	//		Result->Type1_DiscardPacketNum = Result->Type1_DiscardPacketNum + Result->DiscardPacketNum[i];
	//		Result->RateSatisfaction[i] = (((Result->Throughput[i] * 1000000) / simulation_time) / UEList[i].bit_rate) * 100;
	//		if (Result->RateSatisfaction[i] >= 100)
	//			Result->RateSatisfaction[i] = 100;
	//		Result->DelaySatisfaction[i] = ((double)Result->SchedulePackerNum[i] / (double)(Result->DiscardPacketNum[i] + Result->SchedulePackerNum[i])) * 100;
	//	}
	//	Result->Type1_AverageThroughput = Result->Type1_TotalThroughput / DB50_UEnumber;
	//	Result->Type1_AverageDelay = Type1_DelayTemp / DB50_UEnumber;
	//	Result->Type1_PacketLossRatio = ((double)Result->Type1_DiscardPacketNum / (double)(Result->Type1_SchedulePacketNum + Result->Type1_DiscardPacketNum)) * 100;
	//}

	//// 計算type2(Video)的throughput、delay、schedule packet數、discard packet數、rate滿意度、delay滿意度
	//if (DB100_UEnumber > 0)
	//{
	//	for (int i = DB50_UEnumber; i<DB50_UEnumber + DB100_UEnumber; i++)
	//	{
	//		Result->Type2_TotalThroughput = Result->Type2_TotalThroughput + Result->Throughput[i];
	//		Type2_DelayTemp = Type2_DelayTemp + (Result->Delay[i] / Result->SchedulePackerNum[i]);
	//		Result->Type2_SchedulePacketNum = Result->Type2_SchedulePacketNum + Result->SchedulePackerNum[i];
	//		Result->Type2_DiscardPacketNum = Result->Type2_DiscardPacketNum + Result->DiscardPacketNum[i];
	//		Result->RateSatisfaction[i] = (((Result->Throughput[i] * 1000000) / simulation_time) / UEList[i].bit_rate) * 100;
	//		if (Result->RateSatisfaction[i] >= 100)
	//			Result->RateSatisfaction[i] = 100;
	//		Result->DelaySatisfaction[i] = ((double)Result->SchedulePackerNum[i] / (double)(Result->DiscardPacketNum[i] + Result->SchedulePackerNum[i])) * 100;
	//	}
	//	Result->Type2_AverageThroughput = Result->Type2_TotalThroughput / DB100_UEnumber;
	//	Result->Type2_AverageDelay = Type2_DelayTemp / DB100_UEnumber;
	//	Result->Type2_PacketLossRatio = ((double)Result->Type2_DiscardPacketNum / (double)(Result->Type2_SchedulePacketNum + Result->Type2_DiscardPacketNum)) * 100;
	//}

	//// 計算type3的throughput、delay、schedule packet數、discard packet數、rate滿意度、delay滿意度
	//if (DB300_UEnumber > 0)
	//{
	//	for (int i = DB50_UEnumber + DB100_UEnumber; i<DB50_UEnumber + DB100_UEnumber + DB300_UEnumber; i++)
	//	{
	//		Result->Type3_TotalThroughput = Result->Type3_TotalThroughput + Result->Throughput[i];
	//		Type3_DelayTemp = Type3_DelayTemp + (Result->Delay[i] / Result->SchedulePackerNum[i]);
	//		Result->Type3_SchedulePacketNum = Result->Type3_SchedulePacketNum + Result->SchedulePackerNum[i];
	//		Result->Type3_DiscardPacketNum = Result->Type3_DiscardPacketNum + Result->DiscardPacketNum[i];
	//		Result->RateSatisfaction[i] = (((Result->Throughput[i] * 1000000) / simulation_time) / UEList[i].bit_rate) * 100;
	//		if (Result->RateSatisfaction[i] >= 100)
	//			Result->RateSatisfaction[i] = 100;
	//		Result->DelaySatisfaction[i] = ((double)Result->SchedulePackerNum[i] / (double)(Result->DiscardPacketNum[i] + Result->SchedulePackerNum[i])) * 100;
	//	}
	//	Result->Type3_AverageThroughput = Result->Type3_TotalThroughput / DB300_UEnumber;
	//	Result->Type3_AverageDelay = Type3_DelayTemp / DB300_UEnumber;
	//	Result->Type3_PacketLossRatio = ((double)Result->Type3_DiscardPacketNum / (double)(Result->Type3_SchedulePacketNum + Result->Type3_DiscardPacketNum)) * 100;
	//}
	OutputResult("EqualRB", Result);
}

void OutputResult(string Scheme, SimulationResult *Result)
{
	fstream Write_SimulationResultFile;              // 宣告fstream物件，用來存throughput、delay、packet loss ratio
	string SimulationResultFileName;
	fstream Write_RateSatisfactionFile;              // 宣告fstream物件，用來存每個UE的rate滿意度txt
	string RateSatisfactionFileName;
	fstream Write_DelaySatisfactionFile;             // 宣告fstream物件，用來存每個UE的delay滿意度txt
	string DelaySatisfactionFileName;

	int UEID = UEnumber;
	SimulationResultFileName = IntToString(UEID) + "_Simulation Result.txt";
	Write_SimulationResultFile.open(SimulationResultFileName, ios::out | ios::app);
	if (Write_SimulationResultFile.fail())
		cout << "檔案無法開啟" << endl;
	else
	{
		//Write_SimulationResultFile << Scheme << " ";
		Write_SimulationResultFile << (Result->TotalThroughput * 1000 / simulation_time) * 1000 << " " << Result->AverageSystemTime << " " << Result->AvgSystemTime << " " << Result->AvgSystemTime_paper << endl;
		//Write_SimulationResultFile << (Result->AverageThroughput * 1000 / TTI) * 1000 << endl;
		//Write_SimulationResultFile << Result->AverageDelay << endl;
		//Write_SimulationResultFile << Result->PacketLossRatio << endl;
		//		Write_SimulationResultFile << (Result->Type1_TotalThroughput * 1000 / simulation_time) * 1000 << endl;
		//Write_SimulationResultFile << (Result->Type1_AverageThroughput * 1000 / TTI) * 1000 << endl;
		//Write_SimulationResultFile << Result->Type1_AverageDelay << endl;
		//Write_SimulationResultFile << Result->Type1_PacketLossRatio << endl;
		//		Write_SimulationResultFile << (Result->Type2_TotalThroughput * 1000 / simulation_time) * 1000 << endl;
		//Write_SimulationResultFile << (Result->Type2_AverageThroughput * 1000 / TTI) * 1000 << endl;
		//Write_SimulationResultFile << Result->Type2_AverageDelay << endl;
		//Write_SimulationResultFile << Result->Type2_PacketLossRatio << endl;	
		//		Write_SimulationResultFile << (Result->Type3_TotalThroughput * 1000 / simulation_time) * 1000 << endl;
		//Write_SimulationResultFile << (Result->Type3_AverageThroughput * 1000 / TTI) * 1000 << endl;
		//Write_SimulationResultFile << Result->Type3_AverageDelay << endl;
		//Write_SimulationResultFile << Result->Type3_PacketLossRatio << endl;
	}

	//RateSatisfactionFileName = Scheme + "_Rate SatisfactionFile.txt";
	//Write_RateSatisfactionFile.open(RateSatisfactionFileName, ios::out | ios::trunc);
	//if (Write_RateSatisfactionFile.fail())
	//	cout << "檔案無法開啟" << endl;
	//else
	//{
	//	Write_RateSatisfactionFile.setf(ios::fixed, ios::floatfield);
	//	Write_RateSatisfactionFile.precision(4);
	//	for (int i = 0; i<UEnumber; i++)
	//		Write_RateSatisfactionFile << Result->RateSatisfaction[i] << endl;
	//}

	//DelaySatisfactionFileName = Scheme + "_Delay SatisfactionFile.txt";
	//Write_DelaySatisfactionFile.open(DelaySatisfactionFileName, ios::out | ios::trunc);
	//if (Write_DelaySatisfactionFile.fail())
	//	cout << "檔案無法開啟" << endl;
	//else
	//{
	//	Write_DelaySatisfactionFile.setf(ios::fixed, ios::floatfield);
	//	Write_DelaySatisfactionFile.precision(4);
	//	for (int i = 0; i<UEnumber; i++)
	//		Write_DelaySatisfactionFile << Result->DelaySatisfaction[i] << endl;
	//}
}

void Buffer_Status(int t, BufferStatus *Queue, UE *UEList, vector <double> *TempPacketArrivalTime, SimulationResult *Result)
{
	//看看每個UE在這個TTI有沒有資料來
	int TTIPacketCount[UEnumber] = { 0 };				//計算每個UE在每個TTI時來的packet個數
	int TempPacketArrivalTimeID = 0;					//看已經取到這個UE的第幾個PAT了
	for (int i = 0; i < UEnumber; i++)
	{
		TTIPacketCount[i] = 0;
		TempPacketArrivalTimeID = Queue->TempPacketArrivalTimeIndex[i];
		if (!TempPacketArrivalTime[i].empty() && (TempPacketArrivalTimeID < TempPacketArrivalTime[i].size()))
		{
			while (TempPacketArrivalTime[i][TempPacketArrivalTimeID] < t + 1)
			{
				Queue->PacketArrivalTime[i].push_back(TempPacketArrivalTime[i][TempPacketArrivalTimeID]);
				TempPacketArrivalTimeID += 1;
				if (TempPacketArrivalTime[i].empty() || TempPacketArrivalTimeID > TempPacketArrivalTime[i].size() - 1)
					break;
			}
		}
		Result->TotalPacketNum[i] = Result->TotalPacketNum[i] + TTIPacketCount[i];           // 把目前TTI內來的packet數加到總packet數
		Queue->TempPacketArrivalTimeIndex[i] = TempPacketArrivalTimeID;
	}

	// 計算此時TTI每個UE的buffer裡每個packet的HOL delay
	for (int i = 0; i<UEnumber; i++)
	{
		Queue->PacketHOLDelay[i].clear();
		for (int j = 0; j < Queue->PacketArrivalTime[i].size(); j++)
			Queue->PacketHOLDelay[i].push_back((t + 1) - Queue->PacketArrivalTime[i][j]);
	}

	// 處理HOL delay是否超過delay budget，超過的packet要discard掉
	for (int i = 0; i<UEnumber; i++)
	{
		//cout << PacketHOLDelay[i].size() << endl;

		if (!Queue->PacketHOLDelay[i].empty())									//先檢查buffer裡是否有packet要check其HOL delay有無超過budget
			while (Queue->PacketHOLDelay[i][0] > UEList[i].delay_budget)		//packet的HOL delay是否超過delay budget
			{
				//cout << PacketHOLDelay[i][0] << endl;
				if (Queue->HeadPacketSize[i] < UEList[i].packet_size)
					Result->DiscardIncompletePacketNum[i] = Result->DiscardIncompletePacketNum[i] + 1;		//用來計算被砍掉不完整的packet數

				Queue->PacketHOLDelay[i].erase(Queue->PacketHOLDelay[i].begin());							//因為packet的HOL delay超過delay budget，所以要砍掉第一個packet
				Result->DiscardPacketNum[i] = Result->DiscardPacketNum[i] + 1;								//累計discard掉的packet數
				Result->SystemTime[i] = Result->SystemTime[i] + (t + 1) - Queue->PacketArrivalTime[i][0];	///被discard掉的packet在系統的時間
				Queue->PacketArrivalTime[i].erase(Queue->PacketArrivalTime[i].begin());						//也刪掉它在PacketArrivalTime裡記錄的arrival time
				if (Queue->PacketHOLDelay[i].empty())
					break;
			}
		if (!Queue->PacketArrivalTime[i].empty())          // 計算每個UE在packet discard掉後的buffer裡有多少資料量
		{
			Queue->Buffer[i] = (Queue->PacketArrivalTime[i].size() - 1) * UEList[i].packet_size + Queue->HeadPacketSize[i];
			Queue->BeforeScheduleBuffer[i] = Queue->Buffer[i];
		}
		else
		{
			Queue->Buffer[i] = 0.0;
			Queue->BeforeScheduleBuffer[i] = Queue->Buffer[i];
		}
	}
}

void EqualRB(int t, BufferStatus *Queue, UE *UE, SimulationResult *Result)
{
	double InstantRate[UEnumber] = { 0.0 };			// 在t時預計可以拿多少rate
	double Priority = 0.0;							// scheduling時用的priority

	int NumBufferPacket = 0;
	int NumUEBufferHavePacket = 0;
	for (int j = 0; j < UEnumber; j++)
	{
		if (Queue->PacketArrivalTime[j].size() != 0)
			NumUEBufferHavePacket = NumUEBufferHavePacket + 1;
		NumBufferPacket = NumBufferPacket + Queue->PacketArrivalTime[j].size();
	}
		
	if (NumBufferPacket == 0)
		return;

	// 開始幫UE分配傳送時間
	for (int i = 0; i < UEnumber; i++)
	{
		if (Queue->PacketArrivalTime[i].size() == 0)				//如果UE有資料要傳才計算
			continue;

		int CQI = UE[i].CQI;
		double CarryBit = 0.0;
		CarryBit = ap_capacity[UEList[i].CQI - 1] / NumUEBufferHavePacket;			//UE在這個時間內可以傳送多少資料量

		//開始把資料從UE的buffer裡裝進RB裡
		double RBSizeSpace = 0.0;
		int RBAssign = 1;
		RBSizeSpace = CarryBit;
		while (RBAssign)
		{
			if (Queue->HeadPacketSize[i] > RBSizeSpace)	//第一個packet size比RB可攜帶的資料量大
			{
				Queue->HeadPacketSize[i] = Queue->HeadPacketSize[i] - RBSizeSpace;
				//double TransmissionTime = RBSizeSpace / CarryBit;									// Debug用
				Result->SystemTime[i] = Result->SystemTime[i] + RBSizeSpace / CarryBit;
				RBSizeSpace = 0;
				RBAssign = 0;
			}
			else													//第一個packet size比RB可攜帶的資料量小
			{
				RBSizeSpace = RBSizeSpace - Queue->HeadPacketSize[i];
				Result->Delay[i] = Result->Delay[i] + ((t + 1) - Queue->PacketArrivalTime[i][0]);	// 計算每一個packet delay
				//double TransmissionTime = Queue->HeadPacketSize[i] / CarryBit;					// Debug用
				//double WaitingTime = ((t + 1) - Queue->PacketArrivalTime[i][0]);					// Debug用
				Result->SystemTime[i] = Result->SystemTime[i] + ((t + 1) - Queue->PacketArrivalTime[i][0]) + Queue->HeadPacketSize[i] / CarryBit;		// 計算傳送到UE的時間
				Queue->PacketArrivalTime[i].erase(Queue->PacketArrivalTime[i].begin());
				Result->SchedulePackerNum[i] = Result->SchedulePackerNum[i] + 1;
				Queue->PacketHOLDelay[i].erase(Queue->PacketHOLDelay[i].begin());
				if (Queue->PacketArrivalTime[i].empty())
					RBAssign = 0;
				Queue->HeadPacketSize[i] = UE[i].packet_size;
			}
		}
		Result->Throughput[i] = Result->Throughput[i] + ((CarryBit - RBSizeSpace) / 1000000);   // 計算UE的throughput
	}

	// 計算每個UE在這TTI scheduling後的buffer裡有多少資料量
	for (int i = 0; i<UEnumber; i++)
	{
		if (!Queue->PacketArrivalTime[i].empty())
			Queue->Buffer[i] = (Queue->PacketArrivalTime[i].size() - 1) * UE[i].packet_size + Queue->HeadPacketSize[i];
		else
			Queue->Buffer[i] = 0;
	}
}

int main()
{
	for (int times = 0; times < 20; times++)
	{
		for (int i = 0; i < UEnumber; i++)
		{
			//Traffic request initial
			UEList[i].bit_rate = 10;
			UEList[i].packet_size = 800;
			if (i < UEnumber *0.33)
			{
				DB50_UEnumber++;
				UEList[i].delay_budget = 50;
			}
			else
			{
				if (i < UEnumber *0.66)
				{
					DB100_UEnumber++;
					UEList[i].delay_budget = 100;
				}
				else
				{
					DB300_UEnumber++;
					UEList[i].delay_budget = 300;
				}
			}

			//Coordiante initial
			uniformdistribution(&UEList[i]);

			//Other calculation
			UEList[i].CQI = getCQI(&UEList[i]);
			UEList[i].lambdai = UEList[i].bit_rate / UEList[i].packet_size;
		}

		//give packet arrival time
		srand((unsigned)time(NULL));			//亂數種子
		string FileName;						//檔案名稱
		fstream WriteFile;						//宣告fstream物件
		double BufferTimer = 0.0;				//每個UE在eNB裡對應buffer的時間軸
		double InterArrivalTime = 0.0;			//packet的inter-arrival time
		int AcrossTTI = 0;						//用來判斷UE的時間軸，packet的inter-arrival time有無跨過此TTI
		for (int i = 0; i < UEnumber; i++)
		{
			//cout << "UE" << i << endl;
			BufferTimer = 0.0;
			AcrossTTI = 0;
			string UEIndex = IntToString(i);
			FileName = "UE" + UEIndex + "_PAT.txt";				//PAT=packet arrival time
			WriteFile.open(FileName, ios::out | ios::trunc);	//記錄每個UE的PAT
			if (WriteFile.fail())
				cout << "檔案無法開啟" << endl;
			else
			{
				for (int t = 0; t < simulation_time; t++)
				{
					int TTIPacketCount = 0;
					//計算每個packet的到來時間點，並記錄此時TTI每個UE的buffer量
					while (BufferTimer <= t + 1)				//用來計算此TTI來了幾個packet和此TTI結束時目前buffer裡的資料量
					{
						WriteFile.setf(ios::fixed, ios::floatfield);
						WriteFile.precision(3);
						if (AcrossTTI)							//AcrossTTI = 1為inter arrival time有跨過此TTI; AcrossTTI=0為無
						{
							if (outputPAT == 1)
								WriteFile << BufferTimer << endl;	//記錄每個packet的arrival time
																	//TTIPacketCount++;
						}
						else
						{
							InterArrivalTime = exponentially_Distributed(UEList[i].lambdai);//亂數產生inter-arrival time
							BufferTimer = BufferTimer + InterArrivalTime;                   //紀錄每個UE的時間軸
						}
						if (BufferTimer > t + 1)				//BufferTimer有無超過目前此TTI
						{
							AcrossTTI = 1;
							break;
						}
						else
							if (AcrossTTI)
								AcrossTTI = 0;
							else
							{
								if (outputPAT == 1)
									WriteFile << BufferTimer << endl;	// 記錄每個packet的arrival time
																		//TTIPacketCount++;
							}
						//cout << "Packet arrival time：" << BufferTimer << endl;
					}
					//cout << "第" << t+1 << "個TTI的Packet數：" << TTIPacketCount << endl;
				}
			}
			WriteFile.close();
		}
		cout << "Give PAT end." << endl;

		//讀取所有UE的PAT先暫存起來
		string UEPacketPatternFileName;						//UE packet pattern的檔案名稱
		fstream ReadUEPAT;									//宣告fstream物件
		char UEPacketArrivalTime[200];						//用來佔存txt每一行的資料
		double ArrivalTime = 0.0;							//用來佔存抓出來每一行的資料
		vector <double> TempPacketArrivalTime[UEnumber];		//用來暫存UE的packet pattern
		for (int i = 0; i < UEnumber; i++)
			TempPacketArrivalTime[i].clear();
		int NumUETemp = UEnumber;
		string NumUEIndex = IntToString(NumUETemp);
		for (int i = 0; i < UEnumber; i++)
		{
			string UEIndex = IntToString(i);
			UEPacketPatternFileName = "UE" + UEIndex + "_PAT.txt";
			ReadUEPAT.open(UEPacketPatternFileName, ios::in);
			if (!ReadUEPAT)
				cout << "檔案無法開啟" << endl;
			else
			{
				while (ReadUEPAT >> UEPacketArrivalTime)
				{
					ArrivalTime = atof(UEPacketArrivalTime);
					TempPacketArrivalTime[i].push_back(ArrivalTime);
				}
			}
			ReadUEPAT.close();
		}

		//Simulation start
		BufferStatus EqualRB_Buffer;
		SimulationResult EqualRB_Result;
		for (int i = 0; i < UEnumber; i++)
		{
			EqualRB_Buffer.PacketArrivalTime[i].clear();
			EqualRB_Buffer.PacketHOLDelay[i].clear();
			EqualRB_Buffer.HeadPacketSize[i] = UEList[i].packet_size;
		}

		for (int t = 0; t < simulation_time; t++)
		{
			double processing = 0;
			processing = (double)t / (double)simulation_time * 100;

			if (t % (simulation_time / 20) == 0)
				cout << (double)t / (double)simulation_time * 100 << "%" << endl;
			Buffer_Status(t, &EqualRB_Buffer, UEList, TempPacketArrivalTime, &EqualRB_Result);
			EqualRB(t, &EqualRB_Buffer, UEList, &EqualRB_Result);
		}

		Simulation_Result(UEList, &EqualRB_Result);
	}
}