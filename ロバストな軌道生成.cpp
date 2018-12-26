#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define pai 3.14159265359	
#define dt  0.002
#define Te 2.5      //sec
#define φe ( pai/2.0 ) //rad

#define No 100//魚の数
#define T_max 200 //最大反復回数
#define N 4   //次数
//#define r 0.05 //探索範囲
#define w 0.72
#define c1 2.05
#define c2 2.05 //重み

#define L1 0.510
#define L0 0.608
#define m 0.127
#define M 0.035

#define Δζ ( 3.140e-4 )
#define Δη ( -3.470e-4 )
#define c ( 3.206e-3 )
#define k ( 5.349e-3 )

#define ζbar ( (M*L0*L0/12.0) + Δζ )
#define ηbar ( (M*L0*L0/12.0) + Δη )

double φ2_max = 0.0;
int count = 0;
//u(t)
double U0f(double t, double an[N]) {
	int n;
	double sigma_sin;

	sigma_sin = 0.0;
	for (n = 0; n < N; n++) 	sigma_sin += an[n] * (sin(2.0 * (n + 1) * pai * t / Te));

	return (t / Te + sigma_sin);
}

//u(t)ビブン
double U1f(double t, double an[N]) {

	int n;
	double sigma_cos;

	sigma_cos = 0.0;

	for (n = 0; n < N; n++) sigma_cos += an[n] * (n + 1) * (cos(2.0 * (n + 1) * pai * t / Te));

	return ((1 / Te) *(1.0 + 2.0 * pai*sigma_cos));
}

//u(t)ニカイビブン
double U2f(double t, double an[N]) {
	int n;
	double sigma_sin2;

	sigma_sin2 = 0.0;

	for (n = 0; n < N; n++) sigma_sin2 += -(an[n] * (n + 1) * (n + 1) * (sin(2.0 * (n + 1) * pai * t / Te)));

	return (pow((2.0*pai) / Te, 2)*sigma_sin2);
}




// サイクロイド関数 
double φ0f(double U0) {


	return(φe * (U0 - (sin(2.0 * pai * U0) / (2.0 * pai))));

}

double φ1f(double U0, double U1) {

	return(φe * U1 * (1.0 - cos(2.0*pai*U0)));
}

double φ2f(double U0, double U1, double U2) {

	return(φe * (U2 - U2 * cos(2.0*pai*U0) + pow(U1, 2) * 2.0 * pai * sin(2.0*pai*U0)));
}


//θ1についての式
double θ1f1(double θ12) {

	return(θ12);
}

double θ1f2(double θ11, double θ21, double θ22, double φ1, double φ2) {

	double  term1, term2, term3, term4, term5;

	double α1, α2, α3, α4;

	double L = 0.5;

	α1 = ((M*L0*L0 / 3.0) + m*L*L);
	α2 = (((M*L0 / 2.0) + m*L) * L1);
	α3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	α4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = (α4 - Δη + α2 * sin(θ11)) * sin(θ21) * φ2;
	term2 = -(α1 + Δζ) * sin(θ11) * cos(θ11) * θ22 * θ22;
	term3 = -α3 * sin(θ11) * cos(θ21);
	term4 = (2.0 * (α4 - Δζ) * cos(θ11) * cos(θ11) + (Δζ - Δη)) * cos(θ21) * φ1 * θ22;
	term5 = (α2 + (α1 + Δζ) * sin(θ11) * cos(θ21) * cos(θ21)) * cos(θ11) * φ1 * φ1;

	//printf("\tterm1 = %8.4f\tterm2 = %8.4f\tterm3 = %8.4f\tterm4 = %8.4f\tterm5 = %8.4f\n", term1, term2, term3, term4, term5);
	//printf("%10.3e\n", (term1 + term2 + term3 + term4 + term5) / (α1 + Δη));

	return ((term1 + term2 + term3 + term4 + term5) / (α1 + Δη));
}

//θ2についての式
double θ2f1(double θ22) {

	return(θ22);
}

double θ2f2(double θ11, double θ12, double θ21, double θ22, double φ1, double φ2) {

	double term1, term2, term3, term4, term5, term6, term7;

	double α1, α2, α3, α4;
	double L = 0.5;

	α1 = ((M*L0*L0 / 3.0) + m*L*L);
	α2 = (((M*L0 / 2.0) + m*L) * L1);
	α3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	α4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = -(α2 + (α4 - Δζ) * sin(θ11)) * cos(θ11) * cos(θ21) * φ2;
	term2 = -α3 * sin(θ21) * cos(θ11);
	term3 = -c * θ22;
	term4 = -k * θ21;
	term5 = ((α1 + Δζ) * cos(θ11) * cos(θ11) - (Δζ - Δη)) * cos(θ21) * sin(θ21) * φ1 * φ1;
	term6 = -(2.0 * (α4 - Δζ) * cos(θ11) * cos(θ11) + (Δζ - Δη)) * cos(θ21) * φ1 * θ12;
	term7 = 2.0 * (α1 + Δζ) * sin(θ11) * cos(θ21) * θ12 * θ22;

	return ((term1 + term2 + term3 + term4 + term5 + term6 + term7) / ((α1 + Δζ) * cos(θ11) * cos(θ11)));
}

//L=0.0

//θ1についての式
double θ1f1_B(double θ12) {

	return(θ12);
}

double θ1f2_B(double θ11, double θ21, double θ22, double φ1, double φ2) {

	double  term1, term2, term3, term4, term5;
	double α1, α2, α3, α4;
	double L = 0.1;

	α1 = ((M*L0*L0 / 3.0) + m*L*L);
	α2 = (((M*L0 / 2.0) + m*L) * L1);
	α3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	α4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = (α4 - Δη + α2 * sin(θ11)) * sin(θ21) * φ2;
	term2 = -(α1 + Δζ) * sin(θ11) * cos(θ11) * θ22 * θ22;
	term3 = -α3 * sin(θ11) * cos(θ21);
	term4 = (2.0 * (α4 - Δζ) * cos(θ11) * cos(θ11) + (Δζ - Δη)) * cos(θ21) * φ1 * θ22;
	term5 = (α2 + (α1 + Δζ) * sin(θ11) * cos(θ21) * cos(θ21)) * cos(θ11) * φ1 * φ1;

	return ((term1 + term2 + term3 + term4 + term5) / (α1 + Δη));
}

//θ2についての式
double θ2f1_B(double θ22) {

	return(θ22);
}

double θ2f2_B(double θ11, double θ12, double θ21, double θ22, double φ1, double φ2) {

	double term1, term2, term3, term4, term5, term6, term7;
	double α1, α2, α3, α4;
	double L = 0.1;

	α1 = ((M*L0*L0 / 3.0) + m*L*L);
	α2 = (((M*L0 / 2.0) + m*L) * L1);
	α3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	α4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = -(α2 + (α4 - Δζ) * sin(θ11)) * cos(θ11) * cos(θ21) * φ2;
	term2 = -α3 * sin(θ21) * cos(θ11);
	term3 = -c * θ22;
	term4 = -k * θ21;
	term5 = ((α1 + Δζ) * cos(θ11) * cos(θ11) - (Δζ - Δη)) * cos(θ21) * sin(θ21) * φ1 * φ1;
	term6 = -(2.0 * (α4 - Δζ) * cos(θ11) * cos(θ11) + (Δζ - Δη)) * cos(θ21) * φ1 * θ12;
	term7 = 2.0 * (α1 + Δζ) * sin(θ11) * cos(θ21) * θ12 * θ22;

	return ((term1 + term2 + term3 + term4 + term5 + term6 + term7) / ((α1 + Δζ) * cos(θ11) * cos(θ11)));
}



//ルンゲクッタの関数（評価値の計算）
double runge(double an[N], FILE *fp, int print) {
	int  n;
	double sigma_sin, sigma_cos, sigma_sin2;
	double θ11, θ12, θ21, θ22;
	double θ11B, θ12B, θ21B, θ22B;
	double θ11_max, θ21_max, θ11B_max, θ21B_max;
	double φ0, φ1, φ2;
	double t;
	double theta1_1, theta1_2, theta2_1, theta2_2;
	double theta1_1B, theta1_2B, theta2_1B, theta2_2B;
	double d11, d12, d13, d14, d21, d22, d23, d24;
	double d11B, d12B, d13B, d14B, d21B, d22B, d23B, d24B;
	double b11, b12, b13, b14, b21, b22, b23, b24;
	double b11B, b12B, b13B, b14B, b21B, b22B, b23B, b24B;
	double sum = 0.0;
	double U0, U1, U2;

	θ11 = θ12 = θ21 = θ22 = θ11B = θ12B = θ21B = θ22B = 0.0;
	φ0 = φ1 = φ2 = 0.0;
	θ11_max = θ21_max = θ11B_max = θ21B_max = 0.0;

	for (t = 0.0; t <= (Te + 2.0); t += dt)
	{
		if (t <= Te)
		{
			U0 = U0f(t, an);
			U1 = U1f(t, an);
			U2 = U2f(t, an);

			φ0 = φ0f(U0);
			φ1 = φ1f(U0, U1);
			φ2 = φ2f(U0, U1, U2);

			if (fabs(φ2) > 8) {
				return(10000);
			}
			if (print == 1 && φ2 > φ2_max) {
				φ2_max = φ2;
				//printf("%f\n", φ2_max);
			}
					
		}
		else {
			U0 = 1.0;
			U1 = 0.0;
			U2 = 0.0;

			φ0 = φe;
			φ1 = 0;
			φ2 = 0;
		}

		//ルンゲクッタ
		b11 = θ2f1(θ22) * dt;
		b21 = θ2f2(θ11, θ12, θ21, θ22, φ1, φ2) * dt;
		d11 = θ1f1(θ12) * dt;
		d21 = θ1f2(θ11, θ21, θ22, φ1, φ2) * dt;

		theta1_1 = θ11 + d11 * 0.5;
		theta1_2 = θ12 + d21 * 0.5;
		theta2_1 = θ21 + b11 * 0.5;
		theta2_2 = θ22 + b21 * 0.5;

		d12 = θ1f1(theta1_2) * dt;
		d22 = θ1f2(theta1_1, theta2_1, theta2_2, φ1, φ2) * dt;
		b12 = θ2f1(theta2_2) * dt;
		b22 = θ2f2(theta1_1, theta1_2, theta2_1, theta2_2, φ1, φ2) * dt;


		theta1_1 = θ11 + d12 * 0.5;
		theta1_2 = θ12 + d22 * 0.5;
		theta2_1 = θ21 + b12 * 0.5;
		theta2_2 = θ22 + b22 * 0.5;

		d13 = θ1f1(theta1_2) * dt;
		d23 = θ1f2(theta1_1, theta2_1, theta2_2, φ1, φ2) * dt;
		b13 = θ2f1(theta2_2) * dt;
		b23 = θ2f2(theta1_1, theta1_2, theta2_1, theta2_2, φ1, φ2) * dt;

		theta1_1 = θ11 + d13;
		theta1_2 = θ12 + d23;
		theta2_1 = θ21 + b13;
		theta2_2 = θ22 + b23;

		d14 = θ1f1(theta1_2) * dt;
		d24 = θ1f2(theta1_1, theta2_1, theta2_2, φ1, φ2) * dt;
		b14 = θ2f1(theta2_2) * dt;
		b24 = θ2f2(theta1_1, theta1_2, theta2_1, theta2_2, φ1, φ2) * dt;

		θ11 += (d11 + 2.0 * d12 + 2.0 * d13 + d14) / 6.0;
		θ12 += (d21 + 2.0 * d22 + 2.0 * d23 + d24) / 6.0;
		θ21 += (b11 + 2.0 * b12 + 2.0 * b13 + b14) / 6.0;
		θ22 += (b21 + 2.0 * b22 + 2.0 * b23 + b24) / 6.0;



		//L=0.0
		d11B = θ1f1_B(θ12B) * dt;
		d21B = θ1f2_B(θ11B, θ21B, θ22B, φ1, φ2) * dt;
		b11B = θ2f1_B(θ22B) * dt;
		b21B = θ2f2_B(θ11B, θ12B, θ21B, θ22B, φ1, φ2) * dt;


		theta1_1B = θ11B + d11B * 0.5;
		theta1_2B = θ12B + d21B * 0.5;
		theta2_1B = θ21B + b11B * 0.5;
		theta2_2B = θ22B + b21B * 0.5;

		d12B = θ1f1_B(theta1_2B) * dt;
		d22B = θ1f2_B(theta1_1B, theta2_1B, theta2_2B, φ1, φ2) * dt;
		b12B = θ2f1_B(theta2_2B) * dt;
		b22B = θ2f2_B(theta1_1B, theta1_2B, theta2_1B, theta2_2B, φ1, φ2) * dt;


		theta1_1B = θ11B + d12B * 0.5;
		theta1_2B = θ12B + d22B * 0.5;
		theta2_1B = θ21B + b12B * 0.5;
		theta2_2B = θ22B + b22B * 0.5;

		d13B = θ1f1_B(theta1_2B) * dt;
		d23B = θ1f2_B(theta1_1B, theta2_1B, theta2_2B, φ1, φ2) * dt;
		b13B = θ2f1_B(theta2_2B) * dt;
		b23B = θ2f2_B(theta1_1B, theta1_2B, theta2_1B, theta2_2B, φ1, φ2) * dt;

		theta1_1B = θ11B + d13B;
		theta1_2B = θ12B + d23B;
		theta2_1B = θ21B + b13B;
		theta2_2B = θ22B + b23B;


		d14B = θ1f1(theta1_2B) * dt;
		d24B = θ1f2(theta1_1B, theta2_1B, theta2_2B, φ1, φ2) * dt;
		b14B = θ2f1(theta2_2B) * dt;
		b24B = θ2f2(theta1_1B, theta1_2B, theta2_1B, theta2_2B, φ1, φ2) * dt;

		θ11B += (d11B + 2.0 * d12B + 2.0 * d13B + d14B) / 6.0;
		θ12B += (d21B + 2.0 * d22B + 2.0 * d23B + d24B) / 6.0;
		θ21B += (b11B + 2.0 * b12B + 2.0 * b13B + b14B) / 6.0;
		θ22B += (b21B + 2.0 * b22B + 2.0 * b23B + b24B) / 6.0;


		//printf("%10.3e,%10.3e,%10.3e,%10.3e\n", θ11, θ21, θ11B, θ21B);

		if (θ11 > 999 && θ21 > 999) {
			printf("エラー");
			printf("%10.3e,%10.3e\n", θ11, θ21);
		}

		if (t >= Te) {
			if (fabs(θ11) >= θ11_max) θ11_max = fabs(θ11);
			if (fabs(θ21) >= θ21_max)	θ21_max = fabs(θ21);
			if (fabs(θ11B) >= θ11B_max) θ11B_max = fabs(θ11B);
			if (fabs(θ21B) >= θ21B_max) θ21B_max = fabs(θ21B);

			//printf("θ11_max=%10.3e,θ21_max=%10.3e,θ11B_max=%10.3e,θ21_max=%10.3e\n", θ11_max, θ21_max, θ11B_max, θ21B_max);

		}

		if (print == 1)
			fprintf(fp, "%.3f,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n", t, φ0, φ1, φ2, θ11*(180.0 / pai), θ21*(180.0 / pai), θ11B * (180.0 / pai), θ21B * (180.0 / pai));

	}
	
	sum = θ11_max + θ21_max + θ11B_max + θ21B_max;

	if (sum > 99999) {
		puts("エラー!");
		printf("%10.3e\n", sum);
		printf("%10.3e,%10.3e,%10.3e,%10.3e\n", θ11_max, θ21_max, θ11B_max, θ21B_max);

	}

	return (sum);		//評価関数

}

void PSO(FILE *fp) {

	fprintf(fp, ",,,,L=0.5,L=0.1,\n");
	fprintf(fp, "t,φ,φ',φ'',θ1,θ2,θ1,θ2,\n");

	int i, j, co;
	double x[No][N]/*位置ベクトル*/, v[No][N]/*速度ベクトル*/;
	double pbest[No][N + 1]/*最良解と評価値*/;
	double gbest[N + 1]/*全体の最良解*/, fg/*全体の評価値*/;
	double an[N];
	double min = 99999; //初期評価値（出来るだけ大きく）

	srand((unsigned int)time(NULL)); //時間を使った乱数

									 /*Step1 : 初期化*/

	for (i = 0; i < No; i++) {

		for (j = 0; j < N; j++) {

			if (j == 0) {
				/*初期位置と初期速度をランダムに与えます*/
				x[i][j] = (((double)rand() / 36727.0 * 2.0) - 1.0) * -0.1;		//-1.0~1.0
				v[i][j] = ((((double)rand() / 36727.0 * 2.0) - 1.0) * 0.1) / 15.0;

				pbest[i][j] = an[j] = x[i][j];
			}

			else {
				/*初期位置と初期速度をランダムに与えます*/
				x[i][j] = (((double)rand() / 36727.0 * 2.0) - 1.0) * -0.05;		//-1.0~1.0
				v[i][j] = ((((double)rand() / 36727.0 * 2.0) - 1.0) * 0.05) / 15.0;

				pbest[i][j] = an[j] = x[i][j];

			}
		}

		pbest[i][N] = runge(an, fp, 0);				//評価値の計算

		if (pbest[i][N] < min) {			    //もし，評価値が初期評価値(99999)よりも小さければ，
			min = gbest[N] = pbest[i][N];		//minを現在の評価値(pbest[i][N])に置き換えます(min = gbest)
			count++;
			for (j = 0; j < N; j++) gbest[j] = pbest[i][j];
		}
	}
	puts("Step1\n");
	for (j = 0; j < N; j++) {
		printf("an[%d] = %10.3e\n", j, gbest[j]);
		//printf("an[%d] = %f\n", j, gbest[j]);
	}
	printf("評価値 = %10.3e\n\n", min);

	/*Step2 : 位置と速度の更新*/
	for (co = 0; co < T_max; co++) {		//カウンタ
		printf("\r%dかいめ", co + 1);

		for (i = 0; i < No; i++) {
			fg = 0.0;
			for (j = 0; j < N; j++) {

				if (j == 0) {
					v[i][j] = w * v[i][j] + c1 * (double)rand() / 36727.0 * (pbest[i][j] - x[i][j]) + c2 * (double)rand() / 36727.0 * (gbest[j] - x[i][j]);
					x[i][j] += v[i][j];

					if (x[i][j] > 0.1) x[i][j] = 0.1;
					else if (x[i][j] < -0.1) x[i][j] = -0.1;

					an[j] = x[i][j];
				}

				else {

					v[i][j] = w * v[i][j] + c1 * (double)rand() / 36727.0 * (pbest[i][j] - x[i][j]) + c2 * (double)rand() / 36727.0 * (gbest[j] - x[i][j]);
					x[i][j] += v[i][j];

					if (x[i][j] > 0.05) x[i][j] = 0.05;
					else if (x[i][j] < -0.05) x[i][j] = -0.05;

					an[j] = x[i][j];

				}
			}

			fg = runge(an, fp, 0);			//新しい評価値の計算

			/*Step3 : pbestとgbestの更新*/
			if (fg < pbest[i][N]) {
				pbest[i][N] = fg;

				for (j = 0; j < N; j++) {
					pbest[i][j] = x[i][j];

					if (fg < gbest[N]) {
						gbest[N] = fg;
						count++;
						for (j = 0; j < N; j++) gbest[j] = x[i][j];
					}
				}
			}
		}
	}

	//更新します！！
	for (i = 0; i < N; i++) an[i] = gbest[i];

	//for (i = 0; i < N; i++) an[i] = 0.0;
	fg = runge(an, fp, 1);


	puts("\nStep3\n");
	for (j = 0; j < N; j++) {
		printf("an[%d] = %f\n", j, an[j]);
		fprintf(fp, "an[%d] = %f,\n", j, an[j]);
	}
	printf("評価値 = %f %d回更新しました\n\t\n", gbest[N], count);
	fprintf(fp, "評価値=%f,\n", gbest[N]);
	fprintf(fp, "φ2_max=%f,\n", φ2_max);
	fclose(fp);
}

int main(void) {

	puts("***************************************\n");

	printf("1回目\n");
	FILE *fp1;
	char *fname1 = "L=0.5,0.1,Te=2.5_1.csv";
	fopen_s(&fp1, fname1, "w");
	PSO(fp1);
	fclose(fp1);


	puts("***************************************\n");

	printf("2回目\n");
	FILE *fp2;
	char *fname2 = "L=0.5,0.1,Te=2.5_2.csv";
	fopen_s(&fp2, fname2, "w");
	PSO(fp2);
	fclose(fp2);


	puts("***************************************\n");

	printf("3回目\n");
	FILE *fp3;
	char *fname3 = "L=0.5,0.1,Te=2.5_3.csv";
	fopen_s(&fp3, fname3, "w");
	PSO(fp3);
	fclose(fp3);


	puts("***************************************\n");

	printf("4回目\n");
	FILE *fp4;
	char *fname4 = "L=0.5,0.1,Te=2.5_4.csv";
	fopen_s(&fp4, fname4, "w");
	PSO(fp4);
	fclose(fp4);

	puts("***************************************\n");

	printf("5回目\n");
	FILE *fp5;
	char *fname5 = "L=0.5,0.1,Te=2.5_5.csv";
	fopen_s(&fp5, fname5, "w");
	PSO(fp5);
	fclose(fp5);

	puts("***************************************\n");

	printf("6回目\n");
	FILE *fp6;
	char *fname6 = "L=0.5,0.1,Te=2.5_6.csv";
	fopen_s(&fp6, fname6, "w");
	PSO(fp6);
	fclose(fp6);


	puts("***************************************\n");

	printf("7回目\n");
	FILE *fp7;
	char *fname7 = "L=0.5,0.1,Te=2.5_7.csv";
	fopen_s(&fp7, fname7, "w");
	PSO(fp7);
	fclose(fp7);

	puts("***************************************\n");

	printf("8回目\n");
	FILE *fp8;
	char *fname8 = "L=0.5,0.1,Te=2.5_8.csv";
	fopen_s(&fp8, fname8, "w");
	PSO(fp8);
	fclose(fp8);

	puts("***************************************\n");

	printf("9回目\n");
	FILE *fp9;
	char *fname9 = "L=0.5,0.1,Te=2.5_9.csv";
	fopen_s(&fp9, fname9, "w");
	PSO(fp9);
	fclose(fp9);



	puts("***************************************\n");

	printf("10回目\n");
	FILE *fp10;
	char *fname10 = "L=0.5,0.1,Te=2.5_10.csv";
	fopen_s(&fp10, fname10, "w");
	PSO(fp10);
	fclose(fp10);

}