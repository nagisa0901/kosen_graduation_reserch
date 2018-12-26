#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define pai 3.14159265359	
#define dt  0.002
#define Te 2.5      //sec
#define ��e ( pai/2.0 ) //rad

#define No 100//���̐�
#define T_max 200 //�ő唽����
#define N 4   //����
//#define r 0.05 //�T���͈�
#define w 0.72
#define c1 2.05
#define c2 2.05 //�d��

#define L1 0.510
#define L0 0.608
#define m 0.127
#define M 0.035

#define ���� ( 3.140e-4 )
#define ���� ( -3.470e-4 )
#define c ( 3.206e-3 )
#define k ( 5.349e-3 )

#define ��bar ( (M*L0*L0/12.0) + ���� )
#define ��bar ( (M*L0*L0/12.0) + ���� )

double ��2_max = 0.0;
int count = 0;
//u(t)
double U0f(double t, double an[N]) {
	int n;
	double sigma_sin;

	sigma_sin = 0.0;
	for (n = 0; n < N; n++) 	sigma_sin += an[n] * (sin(2.0 * (n + 1) * pai * t / Te));

	return (t / Te + sigma_sin);
}

//u(t)�r�u��
double U1f(double t, double an[N]) {

	int n;
	double sigma_cos;

	sigma_cos = 0.0;

	for (n = 0; n < N; n++) sigma_cos += an[n] * (n + 1) * (cos(2.0 * (n + 1) * pai * t / Te));

	return ((1 / Te) *(1.0 + 2.0 * pai*sigma_cos));
}

//u(t)�j�J�C�r�u��
double U2f(double t, double an[N]) {
	int n;
	double sigma_sin2;

	sigma_sin2 = 0.0;

	for (n = 0; n < N; n++) sigma_sin2 += -(an[n] * (n + 1) * (n + 1) * (sin(2.0 * (n + 1) * pai * t / Te)));

	return (pow((2.0*pai) / Te, 2)*sigma_sin2);
}




// �T�C�N���C�h�֐� 
double ��0f(double U0) {


	return(��e * (U0 - (sin(2.0 * pai * U0) / (2.0 * pai))));

}

double ��1f(double U0, double U1) {

	return(��e * U1 * (1.0 - cos(2.0*pai*U0)));
}

double ��2f(double U0, double U1, double U2) {

	return(��e * (U2 - U2 * cos(2.0*pai*U0) + pow(U1, 2) * 2.0 * pai * sin(2.0*pai*U0)));
}


//��1�ɂ��Ă̎�
double ��1f1(double ��12) {

	return(��12);
}

double ��1f2(double ��11, double ��21, double ��22, double ��1, double ��2) {

	double  term1, term2, term3, term4, term5;

	double ��1, ��2, ��3, ��4;

	double L = 0.5;

	��1 = ((M*L0*L0 / 3.0) + m*L*L);
	��2 = (((M*L0 / 2.0) + m*L) * L1);
	��3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	��4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = (��4 - ���� + ��2 * sin(��11)) * sin(��21) * ��2;
	term2 = -(��1 + ����) * sin(��11) * cos(��11) * ��22 * ��22;
	term3 = -��3 * sin(��11) * cos(��21);
	term4 = (2.0 * (��4 - ����) * cos(��11) * cos(��11) + (���� - ����)) * cos(��21) * ��1 * ��22;
	term5 = (��2 + (��1 + ����) * sin(��11) * cos(��21) * cos(��21)) * cos(��11) * ��1 * ��1;

	//printf("\tterm1 = %8.4f\tterm2 = %8.4f\tterm3 = %8.4f\tterm4 = %8.4f\tterm5 = %8.4f\n", term1, term2, term3, term4, term5);
	//printf("%10.3e\n", (term1 + term2 + term3 + term4 + term5) / (��1 + ����));

	return ((term1 + term2 + term3 + term4 + term5) / (��1 + ����));
}

//��2�ɂ��Ă̎�
double ��2f1(double ��22) {

	return(��22);
}

double ��2f2(double ��11, double ��12, double ��21, double ��22, double ��1, double ��2) {

	double term1, term2, term3, term4, term5, term6, term7;

	double ��1, ��2, ��3, ��4;
	double L = 0.5;

	��1 = ((M*L0*L0 / 3.0) + m*L*L);
	��2 = (((M*L0 / 2.0) + m*L) * L1);
	��3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	��4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = -(��2 + (��4 - ����) * sin(��11)) * cos(��11) * cos(��21) * ��2;
	term2 = -��3 * sin(��21) * cos(��11);
	term3 = -c * ��22;
	term4 = -k * ��21;
	term5 = ((��1 + ����) * cos(��11) * cos(��11) - (���� - ����)) * cos(��21) * sin(��21) * ��1 * ��1;
	term6 = -(2.0 * (��4 - ����) * cos(��11) * cos(��11) + (���� - ����)) * cos(��21) * ��1 * ��12;
	term7 = 2.0 * (��1 + ����) * sin(��11) * cos(��21) * ��12 * ��22;

	return ((term1 + term2 + term3 + term4 + term5 + term6 + term7) / ((��1 + ����) * cos(��11) * cos(��11)));
}

//L=0.0

//��1�ɂ��Ă̎�
double ��1f1_B(double ��12) {

	return(��12);
}

double ��1f2_B(double ��11, double ��21, double ��22, double ��1, double ��2) {

	double  term1, term2, term3, term4, term5;
	double ��1, ��2, ��3, ��4;
	double L = 0.1;

	��1 = ((M*L0*L0 / 3.0) + m*L*L);
	��2 = (((M*L0 / 2.0) + m*L) * L1);
	��3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	��4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = (��4 - ���� + ��2 * sin(��11)) * sin(��21) * ��2;
	term2 = -(��1 + ����) * sin(��11) * cos(��11) * ��22 * ��22;
	term3 = -��3 * sin(��11) * cos(��21);
	term4 = (2.0 * (��4 - ����) * cos(��11) * cos(��11) + (���� - ����)) * cos(��21) * ��1 * ��22;
	term5 = (��2 + (��1 + ����) * sin(��11) * cos(��21) * cos(��21)) * cos(��11) * ��1 * ��1;

	return ((term1 + term2 + term3 + term4 + term5) / (��1 + ����));
}

//��2�ɂ��Ă̎�
double ��2f1_B(double ��22) {

	return(��22);
}

double ��2f2_B(double ��11, double ��12, double ��21, double ��22, double ��1, double ��2) {

	double term1, term2, term3, term4, term5, term6, term7;
	double ��1, ��2, ��3, ��4;
	double L = 0.1;

	��1 = ((M*L0*L0 / 3.0) + m*L*L);
	��2 = (((M*L0 / 2.0) + m*L) * L1);
	��3 = (((M*L0 / 2.0) + m*L) * 9.80665);
	��4 = ((M*L0*L0 / 6.0) + m*L*L);


	term1 = -(��2 + (��4 - ����) * sin(��11)) * cos(��11) * cos(��21) * ��2;
	term2 = -��3 * sin(��21) * cos(��11);
	term3 = -c * ��22;
	term4 = -k * ��21;
	term5 = ((��1 + ����) * cos(��11) * cos(��11) - (���� - ����)) * cos(��21) * sin(��21) * ��1 * ��1;
	term6 = -(2.0 * (��4 - ����) * cos(��11) * cos(��11) + (���� - ����)) * cos(��21) * ��1 * ��12;
	term7 = 2.0 * (��1 + ����) * sin(��11) * cos(��21) * ��12 * ��22;

	return ((term1 + term2 + term3 + term4 + term5 + term6 + term7) / ((��1 + ����) * cos(��11) * cos(��11)));
}



//�����Q�N�b�^�̊֐��i�]���l�̌v�Z�j
double runge(double an[N], FILE *fp, int print) {
	int  n;
	double sigma_sin, sigma_cos, sigma_sin2;
	double ��11, ��12, ��21, ��22;
	double ��11B, ��12B, ��21B, ��22B;
	double ��11_max, ��21_max, ��11B_max, ��21B_max;
	double ��0, ��1, ��2;
	double t;
	double theta1_1, theta1_2, theta2_1, theta2_2;
	double theta1_1B, theta1_2B, theta2_1B, theta2_2B;
	double d11, d12, d13, d14, d21, d22, d23, d24;
	double d11B, d12B, d13B, d14B, d21B, d22B, d23B, d24B;
	double b11, b12, b13, b14, b21, b22, b23, b24;
	double b11B, b12B, b13B, b14B, b21B, b22B, b23B, b24B;
	double sum = 0.0;
	double U0, U1, U2;

	��11 = ��12 = ��21 = ��22 = ��11B = ��12B = ��21B = ��22B = 0.0;
	��0 = ��1 = ��2 = 0.0;
	��11_max = ��21_max = ��11B_max = ��21B_max = 0.0;

	for (t = 0.0; t <= (Te + 2.0); t += dt)
	{
		if (t <= Te)
		{
			U0 = U0f(t, an);
			U1 = U1f(t, an);
			U2 = U2f(t, an);

			��0 = ��0f(U0);
			��1 = ��1f(U0, U1);
			��2 = ��2f(U0, U1, U2);

			if (fabs(��2) > 8) {
				return(10000);
			}
			if (print == 1 && ��2 > ��2_max) {
				��2_max = ��2;
				//printf("%f\n", ��2_max);
			}
					
		}
		else {
			U0 = 1.0;
			U1 = 0.0;
			U2 = 0.0;

			��0 = ��e;
			��1 = 0;
			��2 = 0;
		}

		//�����Q�N�b�^
		b11 = ��2f1(��22) * dt;
		b21 = ��2f2(��11, ��12, ��21, ��22, ��1, ��2) * dt;
		d11 = ��1f1(��12) * dt;
		d21 = ��1f2(��11, ��21, ��22, ��1, ��2) * dt;

		theta1_1 = ��11 + d11 * 0.5;
		theta1_2 = ��12 + d21 * 0.5;
		theta2_1 = ��21 + b11 * 0.5;
		theta2_2 = ��22 + b21 * 0.5;

		d12 = ��1f1(theta1_2) * dt;
		d22 = ��1f2(theta1_1, theta2_1, theta2_2, ��1, ��2) * dt;
		b12 = ��2f1(theta2_2) * dt;
		b22 = ��2f2(theta1_1, theta1_2, theta2_1, theta2_2, ��1, ��2) * dt;


		theta1_1 = ��11 + d12 * 0.5;
		theta1_2 = ��12 + d22 * 0.5;
		theta2_1 = ��21 + b12 * 0.5;
		theta2_2 = ��22 + b22 * 0.5;

		d13 = ��1f1(theta1_2) * dt;
		d23 = ��1f2(theta1_1, theta2_1, theta2_2, ��1, ��2) * dt;
		b13 = ��2f1(theta2_2) * dt;
		b23 = ��2f2(theta1_1, theta1_2, theta2_1, theta2_2, ��1, ��2) * dt;

		theta1_1 = ��11 + d13;
		theta1_2 = ��12 + d23;
		theta2_1 = ��21 + b13;
		theta2_2 = ��22 + b23;

		d14 = ��1f1(theta1_2) * dt;
		d24 = ��1f2(theta1_1, theta2_1, theta2_2, ��1, ��2) * dt;
		b14 = ��2f1(theta2_2) * dt;
		b24 = ��2f2(theta1_1, theta1_2, theta2_1, theta2_2, ��1, ��2) * dt;

		��11 += (d11 + 2.0 * d12 + 2.0 * d13 + d14) / 6.0;
		��12 += (d21 + 2.0 * d22 + 2.0 * d23 + d24) / 6.0;
		��21 += (b11 + 2.0 * b12 + 2.0 * b13 + b14) / 6.0;
		��22 += (b21 + 2.0 * b22 + 2.0 * b23 + b24) / 6.0;



		//L=0.0
		d11B = ��1f1_B(��12B) * dt;
		d21B = ��1f2_B(��11B, ��21B, ��22B, ��1, ��2) * dt;
		b11B = ��2f1_B(��22B) * dt;
		b21B = ��2f2_B(��11B, ��12B, ��21B, ��22B, ��1, ��2) * dt;


		theta1_1B = ��11B + d11B * 0.5;
		theta1_2B = ��12B + d21B * 0.5;
		theta2_1B = ��21B + b11B * 0.5;
		theta2_2B = ��22B + b21B * 0.5;

		d12B = ��1f1_B(theta1_2B) * dt;
		d22B = ��1f2_B(theta1_1B, theta2_1B, theta2_2B, ��1, ��2) * dt;
		b12B = ��2f1_B(theta2_2B) * dt;
		b22B = ��2f2_B(theta1_1B, theta1_2B, theta2_1B, theta2_2B, ��1, ��2) * dt;


		theta1_1B = ��11B + d12B * 0.5;
		theta1_2B = ��12B + d22B * 0.5;
		theta2_1B = ��21B + b12B * 0.5;
		theta2_2B = ��22B + b22B * 0.5;

		d13B = ��1f1_B(theta1_2B) * dt;
		d23B = ��1f2_B(theta1_1B, theta2_1B, theta2_2B, ��1, ��2) * dt;
		b13B = ��2f1_B(theta2_2B) * dt;
		b23B = ��2f2_B(theta1_1B, theta1_2B, theta2_1B, theta2_2B, ��1, ��2) * dt;

		theta1_1B = ��11B + d13B;
		theta1_2B = ��12B + d23B;
		theta2_1B = ��21B + b13B;
		theta2_2B = ��22B + b23B;


		d14B = ��1f1(theta1_2B) * dt;
		d24B = ��1f2(theta1_1B, theta2_1B, theta2_2B, ��1, ��2) * dt;
		b14B = ��2f1(theta2_2B) * dt;
		b24B = ��2f2(theta1_1B, theta1_2B, theta2_1B, theta2_2B, ��1, ��2) * dt;

		��11B += (d11B + 2.0 * d12B + 2.0 * d13B + d14B) / 6.0;
		��12B += (d21B + 2.0 * d22B + 2.0 * d23B + d24B) / 6.0;
		��21B += (b11B + 2.0 * b12B + 2.0 * b13B + b14B) / 6.0;
		��22B += (b21B + 2.0 * b22B + 2.0 * b23B + b24B) / 6.0;


		//printf("%10.3e,%10.3e,%10.3e,%10.3e\n", ��11, ��21, ��11B, ��21B);

		if (��11 > 999 && ��21 > 999) {
			printf("�G���[");
			printf("%10.3e,%10.3e\n", ��11, ��21);
		}

		if (t >= Te) {
			if (fabs(��11) >= ��11_max) ��11_max = fabs(��11);
			if (fabs(��21) >= ��21_max)	��21_max = fabs(��21);
			if (fabs(��11B) >= ��11B_max) ��11B_max = fabs(��11B);
			if (fabs(��21B) >= ��21B_max) ��21B_max = fabs(��21B);

			//printf("��11_max=%10.3e,��21_max=%10.3e,��11B_max=%10.3e,��21_max=%10.3e\n", ��11_max, ��21_max, ��11B_max, ��21B_max);

		}

		if (print == 1)
			fprintf(fp, "%.3f,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n", t, ��0, ��1, ��2, ��11*(180.0 / pai), ��21*(180.0 / pai), ��11B * (180.0 / pai), ��21B * (180.0 / pai));

	}
	
	sum = ��11_max + ��21_max + ��11B_max + ��21B_max;

	if (sum > 99999) {
		puts("�G���[!");
		printf("%10.3e\n", sum);
		printf("%10.3e,%10.3e,%10.3e,%10.3e\n", ��11_max, ��21_max, ��11B_max, ��21B_max);

	}

	return (sum);		//�]���֐�

}

void PSO(FILE *fp) {

	fprintf(fp, ",,,,L=0.5,L=0.1,\n");
	fprintf(fp, "t,��,��',��'',��1,��2,��1,��2,\n");

	int i, j, co;
	double x[No][N]/*�ʒu�x�N�g��*/, v[No][N]/*���x�x�N�g��*/;
	double pbest[No][N + 1]/*�ŗǉ��ƕ]���l*/;
	double gbest[N + 1]/*�S�̂̍ŗǉ�*/, fg/*�S�̂̕]���l*/;
	double an[N];
	double min = 99999; //�����]���l�i�o���邾���傫���j

	srand((unsigned int)time(NULL)); //���Ԃ��g��������

									 /*Step1 : ������*/

	for (i = 0; i < No; i++) {

		for (j = 0; j < N; j++) {

			if (j == 0) {
				/*�����ʒu�Ə������x�������_���ɗ^���܂�*/
				x[i][j] = (((double)rand() / 36727.0 * 2.0) - 1.0) * -0.1;		//-1.0~1.0
				v[i][j] = ((((double)rand() / 36727.0 * 2.0) - 1.0) * 0.1) / 15.0;

				pbest[i][j] = an[j] = x[i][j];
			}

			else {
				/*�����ʒu�Ə������x�������_���ɗ^���܂�*/
				x[i][j] = (((double)rand() / 36727.0 * 2.0) - 1.0) * -0.05;		//-1.0~1.0
				v[i][j] = ((((double)rand() / 36727.0 * 2.0) - 1.0) * 0.05) / 15.0;

				pbest[i][j] = an[j] = x[i][j];

			}
		}

		pbest[i][N] = runge(an, fp, 0);				//�]���l�̌v�Z

		if (pbest[i][N] < min) {			    //�����C�]���l�������]���l(99999)������������΁C
			min = gbest[N] = pbest[i][N];		//min�����݂̕]���l(pbest[i][N])�ɒu�������܂�(min = gbest)
			count++;
			for (j = 0; j < N; j++) gbest[j] = pbest[i][j];
		}
	}
	puts("Step1\n");
	for (j = 0; j < N; j++) {
		printf("an[%d] = %10.3e\n", j, gbest[j]);
		//printf("an[%d] = %f\n", j, gbest[j]);
	}
	printf("�]���l = %10.3e\n\n", min);

	/*Step2 : �ʒu�Ƒ��x�̍X�V*/
	for (co = 0; co < T_max; co++) {		//�J�E���^
		printf("\r%d������", co + 1);

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

			fg = runge(an, fp, 0);			//�V�����]���l�̌v�Z

			/*Step3 : pbest��gbest�̍X�V*/
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

	//�X�V���܂��I�I
	for (i = 0; i < N; i++) an[i] = gbest[i];

	//for (i = 0; i < N; i++) an[i] = 0.0;
	fg = runge(an, fp, 1);


	puts("\nStep3\n");
	for (j = 0; j < N; j++) {
		printf("an[%d] = %f\n", j, an[j]);
		fprintf(fp, "an[%d] = %f,\n", j, an[j]);
	}
	printf("�]���l = %f %d��X�V���܂���\n\t\n", gbest[N], count);
	fprintf(fp, "�]���l=%f,\n", gbest[N]);
	fprintf(fp, "��2_max=%f,\n", ��2_max);
	fclose(fp);
}

int main(void) {

	puts("***************************************\n");

	printf("1���\n");
	FILE *fp1;
	char *fname1 = "L=0.5,0.1,Te=2.5_1.csv";
	fopen_s(&fp1, fname1, "w");
	PSO(fp1);
	fclose(fp1);


	puts("***************************************\n");

	printf("2���\n");
	FILE *fp2;
	char *fname2 = "L=0.5,0.1,Te=2.5_2.csv";
	fopen_s(&fp2, fname2, "w");
	PSO(fp2);
	fclose(fp2);


	puts("***************************************\n");

	printf("3���\n");
	FILE *fp3;
	char *fname3 = "L=0.5,0.1,Te=2.5_3.csv";
	fopen_s(&fp3, fname3, "w");
	PSO(fp3);
	fclose(fp3);


	puts("***************************************\n");

	printf("4���\n");
	FILE *fp4;
	char *fname4 = "L=0.5,0.1,Te=2.5_4.csv";
	fopen_s(&fp4, fname4, "w");
	PSO(fp4);
	fclose(fp4);

	puts("***************************************\n");

	printf("5���\n");
	FILE *fp5;
	char *fname5 = "L=0.5,0.1,Te=2.5_5.csv";
	fopen_s(&fp5, fname5, "w");
	PSO(fp5);
	fclose(fp5);

	puts("***************************************\n");

	printf("6���\n");
	FILE *fp6;
	char *fname6 = "L=0.5,0.1,Te=2.5_6.csv";
	fopen_s(&fp6, fname6, "w");
	PSO(fp6);
	fclose(fp6);


	puts("***************************************\n");

	printf("7���\n");
	FILE *fp7;
	char *fname7 = "L=0.5,0.1,Te=2.5_7.csv";
	fopen_s(&fp7, fname7, "w");
	PSO(fp7);
	fclose(fp7);

	puts("***************************************\n");

	printf("8���\n");
	FILE *fp8;
	char *fname8 = "L=0.5,0.1,Te=2.5_8.csv";
	fopen_s(&fp8, fname8, "w");
	PSO(fp8);
	fclose(fp8);

	puts("***************************************\n");

	printf("9���\n");
	FILE *fp9;
	char *fname9 = "L=0.5,0.1,Te=2.5_9.csv";
	fopen_s(&fp9, fname9, "w");
	PSO(fp9);
	fclose(fp9);



	puts("***************************************\n");

	printf("10���\n");
	FILE *fp10;
	char *fname10 = "L=0.5,0.1,Te=2.5_10.csv";
	fopen_s(&fp10, fname10, "w");
	PSO(fp10);
	fclose(fp10);

}