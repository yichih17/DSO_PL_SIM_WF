#ifndef _DEFIEN_H
#define _DEFIEN_H

//Scenario patameter

#define radius_AP 185			//m

#define power_AP  23.0			//dBm

#define rho_max 0.99			//Max load of BSs

//Algoritm parameter
#define depth_max 2
#define satisfied_TH 75

//Simulation parameter
#define simulation_time 1000000	//ms(TTI
#define UE_dis_mode 1			//0: uniform 1:hotspot
#define UE_type_number 3		//DB = 50, 100, 300ms
#define UEnumber 40

#include<vector>

enum type_bs { macro, ap, ue };
enum type_ue { type1, type2, type3, type4 };
enum type_distribution { uniform, hotspot };

struct BS;
struct UE;
extern std::vector <BS> vbslist;
extern std::vector <UE> vuelist;

int range_ap[];
extern double ap_capacity[8];

/*attribute of BS*/
struct BS
{
	int num;
	type_bs type;
	double coor_X, coor_Y;
	double lambda;
	double systemT;
	double systemT_constraint;
	std::vector <UE*> connectingUE;
};

/*attribute of UE*/
struct UE
{
	int num;
	type_ue type;
	double coor_X, coor_Y;
	double lambdai;
	int delay_budget;
	double packet_size;
	double bit_rate;
	BS* connecting_BS;
	int CQI;
	std::vector <BS*> availBS;
};

struct result
{
	int outage_number;
};

struct connection_status
{
	std::vector <BS> bslist;
	std::vector <UE> uelist;
	int influence;
	int outage_dso;
};

struct hotspot
{
	double coor_x;
	double coor_y;
};

class BufferStatus
{
public:
	int TempPacketArrivalTimeIndex[UEnumber];       // 用來紀錄在TempPacketArrivalTime這個vector抓取第幾個packet了
	std::vector <double> PacketArrivalTime[UEnumber];     // packet的arrival time(ms)
	std::vector <double> PacketHOLDelay[UEnumber];        // packet的HOL delay(ms)
	double Buffer[UEnumber];                        // 每個UE在eNB裡對應buffer的資料量(bits)
	double BeforeScheduleBuffer[UEnumber];          // 排程前buffer的量
	double AfterScheduleBuffer[UEnumber];           // 排程後buffer的量
	double HeadPacketSize[UEnumber];                // 每個UE buffer裡的第一個packet size有多大(bits)
	BufferStatus()                               // 初始化class的變數
	{
		memset(TempPacketArrivalTimeIndex, 0, sizeof(TempPacketArrivalTimeIndex));
		memset(Buffer, 0, sizeof(Buffer));
		memset(BeforeScheduleBuffer, 0, sizeof(BeforeScheduleBuffer));
		memset(AfterScheduleBuffer, 0, sizeof(AfterScheduleBuffer));
		memset(HeadPacketSize, 0, sizeof(HeadPacketSize));
	}
};

class SimulationResult
{
public:
	double PastDataAmount[UEnumber];				// 每個UE在t時過去拿到的資料量
	double AccumulateRate[UEnumber];				// 過去整個指數平均rate歷史          
	int TotalPacketNum[UEnumber];					// 每個UE在eNB裡對應buffer的來的總packete個數
	int DiscardIncompletePacketNum[UEnumber];		// 用來計算被砍掉的不完整的packet有幾個
	int DiscardPacketNum[UEnumber];					// 每個UE在eNB裡對應buffer裡discard的packet數
	int SchedulePackerNum[UEnumber];				// 每個UE在eNB裡對應buffer裡schedule的packet數
													//double PacketDelay[UEnumber];					// 每個UE平均封包delay
	double Throughput[UEnumber];					// 每個UE的throughput
	double Delay[UEnumber];							// 每個UE的封包在buffer裡等待的時間
	double SystemTime[UEnumber];					// 每個UE的封包從進入buffer到傳送完成的時間
	double RateSatisfaction[UEnumber];				// 每個UE的rate satisfaction
	double DelaySatisfaction[UEnumber];				// 每個UE的delay satisfaction
	double AvgSystemTime;						// Queueing算出來的
	double AvgSystemTime_paper;				// Queueing算出來的(原本論文用的算法)
	double TotalThroughput;						// 用來記錄整體系統的throughput
	int TotalSchedulePacketNum;					// 用來記錄整體系統排程多少packet
	int TotalDiscardPacketNum;					// 用來記錄整體系統砍多少packet
	double AverageThroughput;					// 用來記錄所有UE的平均throughput
	double AverageDelay;						// 用來記錄所有UE的平均delay
	double AverageSystemTime;					// 用來記錄所有UE的平均delay
	double PacketLossRatio;						// 用來記錄整體系統的packet loss ratio	
	double Type1_TotalThroughput;				// 用來記錄type1 UE的throughput
	int Type1_SchedulePacketNum;				// 用來記錄type1 UE排程多少packet
	int Type1_DiscardPacketNum;					// 用來記錄type1 UE砍多少packet
	double Type1_AverageThroughput;				// 用來記錄type1 UE的平均throughput
	double Type1_AverageDelay;					// 用來記錄type1 UE的平均delay
	double Type1_PacketLossRatio;				// 用來記錄type1 UE的packet loss ratio
	double Type2_TotalThroughput;				// 用來記錄type2 UE的throughput
	int Type2_SchedulePacketNum;				// 用來記錄type2 UE排程多少packet
	int Type2_DiscardPacketNum;					// 用來記錄type2 UE砍多少packet	
	double Type2_AverageThroughput;				// 用來記錄type2 UE的平均throughput
	double Type2_AverageDelay;					// 用來記錄type2 UE的平均delay
	double Type2_PacketLossRatio;				// 用來記錄type2 UE的packet loss ratio
	double Type3_TotalThroughput;				// 用來記錄type3 UE的throughput
	int Type3_SchedulePacketNum;				// 用來記錄type3 UE排程多少packet
	int Type3_DiscardPacketNum;					// 用來記錄type3 UE砍多少packet	
	double Type3_AverageThroughput;				// 用來記錄type3 UE的平均throughput
	double Type3_AverageDelay;					// 用來記錄type3 UE的平均delay
	double Type3_PacketLossRatio;				// 用來記錄type3 UE的packet loss ratio
	double Type4_TotalThroughput;				// 用來記錄type4 UE的throughput
	int Type4_SchedulePacketNum;				// 用來記錄type4 UE排程多少packet
	int Type4_DiscardPacketNum;					// 用來記錄type4 UE砍多少packet	
	double Type4_AverageThroughput;				// 用來記錄type4 UE的平均throughput
	double Type4_AverageDelay;					// 用來記錄type4 UE的平均delay
	double Type4_PacketLossRatio;				// 用來記錄type4 UE的packet loss ratio
	SimulationResult()						// 初始化class的變數
	{
		memset(PastDataAmount, 0, sizeof(PastDataAmount));
		memset(AccumulateRate, 1, sizeof(AccumulateRate));
		memset(TotalPacketNum, 0, sizeof(TotalPacketNum));
		memset(DiscardIncompletePacketNum, 0, sizeof(DiscardIncompletePacketNum));
		memset(DiscardPacketNum, 0, sizeof(DiscardPacketNum));
		memset(SchedulePackerNum, 0, sizeof(SchedulePackerNum));
		memset(Throughput, 0, sizeof(Throughput));
		memset(Delay, 0, sizeof(Delay));
		memset(SystemTime, 0, sizeof(SystemTime));
		memset(RateSatisfaction, 0, sizeof(RateSatisfaction));
		memset(DelaySatisfaction, 0, sizeof(DelaySatisfaction));
		AvgSystemTime = 0.0;
		AvgSystemTime_paper = 0.0;
		TotalThroughput = 0.0;
		TotalSchedulePacketNum = 0;
		TotalDiscardPacketNum = 0;
		AverageThroughput = 0.0;
		AverageDelay = 0.0;
		AverageSystemTime = 0.0;
		PacketLossRatio = 0.0;
		Type1_TotalThroughput = 0.0;
		Type1_SchedulePacketNum = 0;
		Type1_DiscardPacketNum = 0;
		Type1_AverageThroughput = 0.0;
		Type1_AverageDelay = 0.0;
		Type1_PacketLossRatio = 0.0;
		Type2_TotalThroughput = 0.0;
		Type2_SchedulePacketNum = 0;
		Type2_DiscardPacketNum = 0;
		Type2_AverageThroughput = 0.0;
		Type2_AverageDelay = 0.0;
		Type2_PacketLossRatio = 0.0;
		Type3_TotalThroughput = 0.0;
		Type3_SchedulePacketNum = 0;
		Type3_DiscardPacketNum = 0;
		Type3_AverageThroughput = 0.0;
		Type3_AverageDelay = 0.0;
		Type3_PacketLossRatio = 0.0;
		Type4_TotalThroughput = 0.0;
		Type4_SchedulePacketNum = 0;
		Type4_DiscardPacketNum = 0;
		Type4_AverageThroughput = 0.0;
		Type4_AverageDelay = 0.0;
		Type4_PacketLossRatio = 0.0;
	}
};

#endif // !_DEFIEN_H