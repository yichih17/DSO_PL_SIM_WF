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

// int��string
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

	// �p����骺throughput�Bdelay�Bschedule packet�ơBdiscard packet��
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

	//// �p��typ1(VoIP)��throughput�Bdelay�Bschedule packet�ơBdiscard packet�ơBrate���N�סBdelay���N��
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

	//// �p��type2(Video)��throughput�Bdelay�Bschedule packet�ơBdiscard packet�ơBrate���N�סBdelay���N��
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

	//// �p��type3��throughput�Bdelay�Bschedule packet�ơBdiscard packet�ơBrate���N�סBdelay���N��
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
	fstream Write_SimulationResultFile;              // �ŧifstream����A�ΨӦsthroughput�Bdelay�Bpacket loss ratio
	string SimulationResultFileName;
	fstream Write_RateSatisfactionFile;              // �ŧifstream����A�ΨӦs�C��UE��rate���N��txt
	string RateSatisfactionFileName;
	fstream Write_DelaySatisfactionFile;             // �ŧifstream����A�ΨӦs�C��UE��delay���N��txt
	string DelaySatisfactionFileName;

	int UEID = UEnumber;
	SimulationResultFileName = IntToString(UEID) + "_Simulation Result.txt";
	Write_SimulationResultFile.open(SimulationResultFileName, ios::out | ios::app);
	if (Write_SimulationResultFile.fail())
		cout << "�ɮ׵L�k�}��" << endl;
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
	//	cout << "�ɮ׵L�k�}��" << endl;
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
	//	cout << "�ɮ׵L�k�}��" << endl;
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
	//�ݬݨC��UE�b�o��TTI���S����ƨ�
	int TTIPacketCount[UEnumber] = { 0 };				//�p��C��UE�b�C��TTI�ɨӪ�packet�Ӽ�
	int TempPacketArrivalTimeID = 0;					//�ݤw�g����o��UE���ĴX��PAT�F
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
		Result->TotalPacketNum[i] = Result->TotalPacketNum[i] + TTIPacketCount[i];           // ��ثeTTI���Ӫ�packet�ƥ[���`packet��
		Queue->TempPacketArrivalTimeIndex[i] = TempPacketArrivalTimeID;
	}

	// �p�⦹��TTI�C��UE��buffer�̨C��packet��HOL delay
	for (int i = 0; i<UEnumber; i++)
	{
		Queue->PacketHOLDelay[i].clear();
		for (int j = 0; j < Queue->PacketArrivalTime[i].size(); j++)
			Queue->PacketHOLDelay[i].push_back((t + 1) - Queue->PacketArrivalTime[i][j]);
	}

	// �B�zHOL delay�O�_�W�Ldelay budget�A�W�L��packet�ndiscard��
	for (int i = 0; i<UEnumber; i++)
	{
		//cout << PacketHOLDelay[i].size() << endl;

		if (!Queue->PacketHOLDelay[i].empty())									//���ˬdbuffer�̬O�_��packet�ncheck��HOL delay���L�W�Lbudget
			while (Queue->PacketHOLDelay[i][0] > UEList[i].delay_budget)		//packet��HOL delay�O�_�W�Ldelay budget
			{
				//cout << PacketHOLDelay[i][0] << endl;
				if (Queue->HeadPacketSize[i] < UEList[i].packet_size)
					Result->DiscardIncompletePacketNum[i] = Result->DiscardIncompletePacketNum[i] + 1;		//�Ψӭp��Q�屼�����㪺packet��

				Queue->PacketHOLDelay[i].erase(Queue->PacketHOLDelay[i].begin());							//�]��packet��HOL delay�W�Ldelay budget�A�ҥH�n�屼�Ĥ@��packet
				Result->DiscardPacketNum[i] = Result->DiscardPacketNum[i] + 1;								//�֭pdiscard����packet��
				Result->SystemTime[i] = Result->SystemTime[i] + (t + 1) - Queue->PacketArrivalTime[i][0];	///�Qdiscard����packet�b�t�Ϊ��ɶ�
				Queue->PacketArrivalTime[i].erase(Queue->PacketArrivalTime[i].begin());						//�]�R�����bPacketArrivalTime�̰O����arrival time
				if (Queue->PacketHOLDelay[i].empty())
					break;
			}
		if (!Queue->PacketArrivalTime[i].empty())          // �p��C��UE�bpacket discard���᪺buffer�̦��h�ָ�ƶq
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
	double InstantRate[UEnumber] = { 0.0 };			// �bt�ɹw�p�i�H���h��rate
	double Priority = 0.0;							// scheduling�ɥΪ�priority

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

	// �}�l��UE���t�ǰe�ɶ�
	for (int i = 0; i < UEnumber; i++)
	{
		if (Queue->PacketArrivalTime[i].size() == 0)				//�p�GUE����ƭn�Ǥ~�p��
			continue;

		int CQI = UE[i].CQI;
		double CarryBit = 0.0;
		CarryBit = ap_capacity[UEList[i].CQI - 1] / NumUEBufferHavePacket;			//UE�b�o�Ӯɶ����i�H�ǰe�h�ָ�ƶq

		//�}�l���ƱqUE��buffer�̸˶iRB��
		double RBSizeSpace = 0.0;
		int RBAssign = 1;
		RBSizeSpace = CarryBit;
		while (RBAssign)
		{
			if (Queue->HeadPacketSize[i] > RBSizeSpace)	//�Ĥ@��packet size��RB�i��a����ƶq�j
			{
				Queue->HeadPacketSize[i] = Queue->HeadPacketSize[i] - RBSizeSpace;
				//double TransmissionTime = RBSizeSpace / CarryBit;									// Debug��
				Result->SystemTime[i] = Result->SystemTime[i] + RBSizeSpace / CarryBit;
				RBSizeSpace = 0;
				RBAssign = 0;
			}
			else													//�Ĥ@��packet size��RB�i��a����ƶq�p
			{
				RBSizeSpace = RBSizeSpace - Queue->HeadPacketSize[i];
				Result->Delay[i] = Result->Delay[i] + ((t + 1) - Queue->PacketArrivalTime[i][0]);	// �p��C�@��packet delay
				//double TransmissionTime = Queue->HeadPacketSize[i] / CarryBit;					// Debug��
				//double WaitingTime = ((t + 1) - Queue->PacketArrivalTime[i][0]);					// Debug��
				Result->SystemTime[i] = Result->SystemTime[i] + ((t + 1) - Queue->PacketArrivalTime[i][0]) + Queue->HeadPacketSize[i] / CarryBit;		// �p��ǰe��UE���ɶ�
				Queue->PacketArrivalTime[i].erase(Queue->PacketArrivalTime[i].begin());
				Result->SchedulePackerNum[i] = Result->SchedulePackerNum[i] + 1;
				Queue->PacketHOLDelay[i].erase(Queue->PacketHOLDelay[i].begin());
				if (Queue->PacketArrivalTime[i].empty())
					RBAssign = 0;
				Queue->HeadPacketSize[i] = UE[i].packet_size;
			}
		}
		Result->Throughput[i] = Result->Throughput[i] + ((CarryBit - RBSizeSpace) / 1000000);   // �p��UE��throughput
	}

	// �p��C��UE�b�oTTI scheduling�᪺buffer�̦��h�ָ�ƶq
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
		srand((unsigned)time(NULL));			//�üƺؤl
		string FileName;						//�ɮצW��
		fstream WriteFile;						//�ŧifstream����
		double BufferTimer = 0.0;				//�C��UE�beNB�̹���buffer���ɶ��b
		double InterArrivalTime = 0.0;			//packet��inter-arrival time
		int AcrossTTI = 0;						//�ΨӧP�_UE���ɶ��b�Apacket��inter-arrival time���L��L��TTI
		for (int i = 0; i < UEnumber; i++)
		{
			//cout << "UE" << i << endl;
			BufferTimer = 0.0;
			AcrossTTI = 0;
			string UEIndex = IntToString(i);
			FileName = "UE" + UEIndex + "_PAT.txt";				//PAT=packet arrival time
			WriteFile.open(FileName, ios::out | ios::trunc);	//�O���C��UE��PAT
			if (WriteFile.fail())
				cout << "�ɮ׵L�k�}��" << endl;
			else
			{
				for (int t = 0; t < simulation_time; t++)
				{
					int TTIPacketCount = 0;
					//�p��C��packet����Ӯɶ��I�A�ðO������TTI�C��UE��buffer�q
					while (BufferTimer <= t + 1)				//�Ψӭp�⦹TTI�ӤF�X��packet�M��TTI�����ɥثebuffer�̪���ƶq
					{
						WriteFile.setf(ios::fixed, ios::floatfield);
						WriteFile.precision(3);
						if (AcrossTTI)							//AcrossTTI = 1��inter arrival time����L��TTI; AcrossTTI=0���L
						{
							if (outputPAT == 1)
								WriteFile << BufferTimer << endl;	//�O���C��packet��arrival time
																	//TTIPacketCount++;
						}
						else
						{
							InterArrivalTime = exponentially_Distributed(UEList[i].lambdai);//�üƲ���inter-arrival time
							BufferTimer = BufferTimer + InterArrivalTime;                   //�����C��UE���ɶ��b
						}
						if (BufferTimer > t + 1)				//BufferTimer���L�W�L�ثe��TTI
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
									WriteFile << BufferTimer << endl;	// �O���C��packet��arrival time
																		//TTIPacketCount++;
							}
						//cout << "Packet arrival time�G" << BufferTimer << endl;
					}
					//cout << "��" << t+1 << "��TTI��Packet�ơG" << TTIPacketCount << endl;
				}
			}
			WriteFile.close();
		}
		cout << "Give PAT end." << endl;

		//Ū���Ҧ�UE��PAT���Ȧs�_��
		string UEPacketPatternFileName;						//UE packet pattern���ɮצW��
		fstream ReadUEPAT;									//�ŧifstream����
		char UEPacketArrivalTime[200];						//�ΨӦ��stxt�C�@�檺���
		double ArrivalTime = 0.0;							//�ΨӦ��s��X�ӨC�@�檺���
		vector <double> TempPacketArrivalTime[UEnumber];		//�ΨӼȦsUE��packet pattern
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
				cout << "�ɮ׵L�k�}��" << endl;
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