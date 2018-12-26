#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<Windows.h>
#include<time.h>

void tansaku(FILE *fp, FILE *fp_É∆1, FILE *fp_É∆2) {

	float É∆1, É∆2;
	float É∆1_max, É∆2_max;
	float t;
	float t_print1, t_print2;
	float iranai1, iranai2, iranai3, iranai4;
	int ret;

	iranai1 = iranai2 = iranai3 = iranai4 = 0.0;
	É∆1_max = É∆2_max = 0.0;
	É∆1 = É∆2 = 0.0;

	while ((ret = fscanf_s(fp, "%f,%f,%f,%f,%f,%f,%f", &t, &iranai1, &iranai2, &É∆2, &É∆1, &iranai3, &iranai4)) != EOF) {

		//printf("%f,\n",t);
		if (t >= 2.5) {
			if (fabs(É∆1) > É∆1_max) {
				É∆1_max = É∆1;
				t_print1 = t;
			}
			if (fabs(É∆2) > É∆2_max) {
				É∆2_max = É∆2;
				t_print2 = t;
			}
		}
	}

	fprintf(fp_É∆1,"%f,%f,\n", t_print1,É∆1_max);
	fprintf(fp_É∆2,"%f,%f,\n", t_print2,É∆2_max);

	printf("%f,%f\n", É∆1_max, É∆2_max);

	printf("èIóπ");
}


int main(void) {

	FILE *fp_É∆1;
	char *fname_É∆1 = "É∆1max.csv";
	fopen_s(&fp_É∆1, fname_É∆1, "w");

	FILE *fp_É∆2;
	char *fname_É∆2 = "É∆2max.csv";
	fopen_s(&fp_É∆2, fname_É∆2, "w");


	FILE *fp1;
	char *fname1 = "optl01_l00.csv";
	fopen_s(&fp1, fname1, "r");
	tansaku(fp1,fp_É∆1,fp_É∆2);
	fclose(fp1);

	FILE *fp2;
	char *fname2 = "optl01_l010.csv";
	fopen_s(&fp2, fname2, "r");
	tansaku(fp2, fp_É∆1, fp_É∆2);
	fclose(fp2);

	FILE *fp3;
	char *fname3 = "optl01_l015.csv";
	fopen_s(&fp3, fname3, "r");
	tansaku(fp3, fp_É∆1, fp_É∆2);
	fclose(fp3);

	FILE *fp4;
	char *fname4 = "optl01_l020.csv";
	fopen_s(&fp4, fname4, "r");
	tansaku(fp4, fp_É∆1, fp_É∆2);
	fclose(fp4);

	FILE *fp5;
	char *fname5 = "optl01_l025.csv";
	fopen_s(&fp5, fname5, "r");
	tansaku(fp5, fp_É∆1, fp_É∆2);
	fclose(fp5);

	FILE *fp6;
	char *fname6 = "optl01_l030.csv";
	fopen_s(&fp6, fname6, "r");
	tansaku(fp6, fp_É∆1, fp_É∆2);
	fclose(fp6);

	FILE *fp7;
	char *fname7 = "optl01_l035.csv";
	fopen_s(&fp7, fname7, "r");
	tansaku(fp7, fp_É∆1, fp_É∆2);
	fclose(fp7);

	FILE *fp8;
	char *fname8 = "optl01_l040.csv";
	fopen_s(&fp8, fname8, "r");
	tansaku(fp8, fp_É∆1, fp_É∆2);
	fclose(fp8);

	FILE *fp9;
	char *fname9 = "optl01_l045.csv";
	fopen_s(&fp9, fname9, "r");
	tansaku(fp9, fp_É∆1, fp_É∆2);
	fclose(fp9);

	FILE *fp10;
	char *fname10 = "optl01_l050.csv";
	fopen_s(&fp10, fname10, "r");
	tansaku(fp10, fp_É∆1, fp_É∆2);
	fclose(fp10);

	fclose(fp_É∆1);
	fclose(fp_É∆2);
}