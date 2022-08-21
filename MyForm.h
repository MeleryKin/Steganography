#pragma comment(lib,"Winmm.lib")
#pragma once
#include <windows.h> 
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <stdio.h>
#include <tchar.h>
#include <conio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

bool nowplaying[5];

namespace Project {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace std;

	bool kontchecked = false, msgchecked = false;
	char kontpath[256], msgpath[256];
	long kontdata, kontsize, msgsize;
	const long n = 900000;
	const long m = n / 16;
	unsigned long byteRate;
	unsigned short audioFormat;

	struct WAVHEADER
	{
		long chunkId;//RIFF
		unsigned long chunkSize;//RIFF size
		long format;//Wave
	};
	WAVHEADER header;
	/// <summary>
	/// Сводка для MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: добавьте код конструктора
			//
		}

		double StandDev(char a[256], string s, double &max) {
			kontsize = kontsize / 2;
			const long n1 = 500000;
			short m1[n1/2], m2[n1/2];
			char b[256];
			for (int i = 0; i < s.length(); i++) {
				b[i] = s[i];
			}
			for (int i = s.length(); i < 256; i++) {
				b[i] = NULL;
			}
			FILE* file;
			errno_t err;
			err = fopen_s(&file, a, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return 0;
			}
			FILE* filenew;
			err = fopen_s(&filenew, b, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return 0;
			}
			fseek(file, kontdata, SEEK_SET);
			fseek(filenew, kontdata, SEEK_SET);
			long long rez;
			max = 0;
			rez = 0;
			for (int i = 0; i < kontsize / n1; i++) {
				fread(&m1, n1, 1, file);
				fread(&m2, n1, 1, filenew);
				for (int j = 0; j < (n1/2); j++) {
					//rez += pow(long(m1[j] - m2[j]), 2);
					long p = (long)m1[j] - (long)m2[j];
					rez += p*p;
					/*if (abs(m2[j]) > max) {
						max = abs(m2[j]);
					}*/
				}
			}
			if (kontsize % n1 > 0) {
				fread(&m1, kontsize % n1, 1, file);
				fread(&m2, kontsize % n1, 1, filenew);
				for (int j = 0; j < (kontsize % n1)/2; j++) {
				//	rez += pow(long(m1[j] - m2[j]), 2);
					long p = (long)m1[j] - (long)m2[j];
					rez += p * p;
					/*if (abs(m2[j]) > max) {
						max = abs(m2[j]);
					}*/
				}
			}
			max = 32767;
		//	MessageBox::Show(Convert::ToString(max));
			kontsize = kontsize * 2;
			fclose(file);
			fclose(filenew);
			double r = (double)rez / kontsize;
			return double (sqrt(r));
		}

		/*double RSAnalysis(string s) {
			double  rs;
			long spl, splflip;
			long RpMp = 0, RmMp = 0, RpMm = 0, RmMm = 0, SpMp = 0, SmMp = 0, SpMm = 0, SmMm = 0;
			short m2[n / 2];
			char b[256];
			for (int i = 0; i < s.length(); i++) {
				b[i] = s[i];
			}
			for (int i = s.length(); i < 256; i++) {
				b[i] = NULL;
			}
			errno_t err;
			FILE* filenew;
			err = fopen_s(&filenew, b, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return 0;
			}
			fseek(filenew, kontdata, SEEK_SET);
			for(int i=0;i<(kontsize/(n/2))/2;i++){
				fread(&m2, n, 1, filenew);
			for (int j = 0; j + 7 < (n / 2); j += 8) {
				short mas[4];
				for (int k = 0; k < 4; k++) {
					mas[k] = m2[j + k * 2];
				}
				spl = 0;
				splflip = 0;
				for (int k = 0; k < 3; k++) {
					spl += abs(mas[k] - mas[k + 1]);
				}
				for (int k = 1; k < 3; k++) {
					if (abs(mas[k]) % 2 == 0) mas[k]++;
					else mas[k]--;
				}
				splflip = 0;
				for (int k = 0; k < 3; k++) {
					splflip += abs(mas[k] - mas[k + 1]);
				}
				if (spl < splflip) SpMp++;
				if (spl > splflip) RpMp++;
				for (int k = 0; k < 4; k++) {
					mas[k] = m2[j + k * 2];
				}
				for (int k = 1; k < 3; k++) {
					if (abs(mas[k]) % 2 == 0) mas[k]--;
					else mas[k]++;
				}
				splflip = 0;
				for (int k = 0; k < 3; k++) {
					splflip += abs(mas[k] - mas[k + 1]);
				}
				if (spl < splflip) SmMp++;
				if (spl > splflip) RmMp++;
				for (int k = 0; k < 4; k++) {
					mas[k] = m2[j + k * 2];
					mas[k] ^= (1 << 0);
				}
				spl = 0;
				for (int k = 0; k < 3; k++) {
					spl += abs(mas[k] - mas[k + 1]);
				}
				for (int k = 1; k < 3; k++) {
					if (abs(mas[k]) % 2 == 0) mas[k]++;
					else mas[k]--;
				}
				splflip = 0;
				for (int k = 0; k < 3; k++) {
					splflip += abs(mas[k] - mas[k + 1]);
				}
				if (spl < splflip) SpMm++;
				if (spl > splflip) RpMm++;
				for (int k = 0; k < 4; k++) {
					mas[k] = m2[j + k * 2];
					mas[k] ^= (1 << 0);
				}
				for (int k = 1; k < 3; k++) {
					if (abs(mas[k]) % 2 == 0) mas[k]--;
					else mas[k]++;
				}
				splflip = 0;
				for (int k = 0; k < 3; k++) {
					splflip += abs(mas[k] - mas[k + 1]);
				}
				if (spl < splflip) SmMm++;
				if (spl > splflip) RmMm++;
			}
		}
			if ((kontsize / (n / 2)) / 2==0) {
				fread(&m2, kontsize*2, 1, filenew);
				for (int j = 0; j + 7 < kontsize; j += 8) {
					short mas[4];
					for (int k = 0; k < 4; k++) {
						mas[k] = m2[j + k * 2];
					}
					spl = 0;
					for (int k = 0; k < 3; k++) {
						spl += abs(mas[k] - mas[k + 1]);
					}
					for (int k = 1; k < 3; k++) {
						if (abs(mas[k]) % 2 == 0) mas[k]++;
						else mas[k]--;
					}
					splflip = 0;
					for (int k = 0; k < 3; k++) {
						splflip += abs(mas[k] - mas[k + 1]);
					}
					if (spl < splflip) SpMp++;
					if (spl > splflip) RpMp++;
					for (int k = 0; k < 4; k++) {
						mas[k] = m2[j + k * 2];
					}
					for (int k = 1; k < 3; k++) {
						if (abs(mas[k]) % 2 == 0) mas[k]--;
						else mas[k]++;
					}
					splflip = 0;
					for (int k = 0; k < 3; k++) {
						splflip += abs(mas[k] - mas[k + 1]);
					}
					if (spl < splflip) SmMp++;
					if (spl > splflip) RmMp++;
					for (int k = 0; k < 4; k++) {
						mas[k] = m2[j + k * 2];
						mas[k] ^= (1 << 0);
					}
					spl = 0;
					for (int k = 0; k < 3; k++) {
						spl += abs(mas[k] - mas[k + 1]);
					}
					for (int k = 1; k < 3; k++) {
						if (abs(mas[k]) % 2 == 0) mas[k]++;
						else mas[k]--;
					}
					splflip = 0;
					for (int k = 0; k < 3; k++) {
						splflip += abs(mas[k] - mas[k + 1]);
					}
					if (spl < splflip) SpMm++;
					if (spl > splflip) RpMm++;
					for (int k = 0; k < 4; k++) {
						mas[k] = m2[j + k * 2];
						mas[k] ^= (1 << 0);
					}
					for (int k = 1; k < 3; k++) {
						if (abs(mas[k]) % 2 == 0) mas[k]--;
						else mas[k]++;
					}
					splflip = 0;
					for (int k = 0; k < 3; k++) {
						splflip += abs(mas[k] - mas[k + 1]);
					}
					if (spl < splflip) SmMm++;
					if (spl > splflip) RmMm++;
				}
			}
			long dp0, dm0, dp1, dm1;
			dp0 = abs(RpMp - SpMp);
			dm0 = abs(RmMp - SmMp);
			dp1 = abs(RpMm - SpMm);
			dm1 = abs(RmMm - SmMm);
			double dis = dm0 - dm1 - dp1 - 3 * dp0;
			double x1, x2;
			if (dp1 + dp0 == 0) {
				x1 = double(-1 * (double(dp0 - dm0)) / dis);
			}
			else {
				dis = dis * dis-8*(dp1 + dp0)*(dp0 - dm0);
				x1 = (double)((-1 * (dm0 - dm1 - dp1 - 3 * dp0) + sqrt(double(dis))) / (4 * (dp1 + dp0)));
				x2 = (double)((-1 * (dm0 - dm1 - dp1 - 3 * dp0) - sqrt(double(dis))) / (4 * (dp1 + dp0)));
				if (abs(x1) > abs(x2))x1 = x2;
			}
			rs = x1 / (x1 - 0.5);
			kontsize = kontsize * 2;
			fclose(filenew);
			return abs(rs);
		}*/

		long NumUnits(long n) {
			long count = 0;
			for (; n; count++)
				n &= (n - 1);
			return count;
		}

		void FileCopy() 
		{
			unsigned char a[n];
			unsigned long d = header.chunkSize;
			FILE* file;
			errno_t err;
			err = fopen_s(&file, kontpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			remove("orig.wav");
			FILE* filenew;
			err = fopen_s(&filenew, "orig.wav", "wb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 1; i <= (d / n); i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, d % n, 1, file);
			fwrite(&a, d % n, 1, filenew);
			fclose(file);
			fclose(filenew);
		}

		void AddEnd() 
		{
			unsigned char a[n];
			unsigned long d = header.chunkSize;
			FILE* file;
			errno_t err;
			err = fopen_s(&file, kontpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			remove("addend.wav");
			FILE* filenew;
			err = fopen_s(&filenew, "addend.wav", "wb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 1; i <= (d / n); i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, d % n, 1, file);
			fwrite(&a, d % n, 1, filenew);
			fclose(file);
			err = fopen_s(&file, msgpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 1; i <= (msgsize / n); i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, msgsize % n, 1, file);
			fwrite(&a, msgsize % n, 1, filenew);
			fclose(file);
			fclose(filenew);
		}

		void LSB() {
			unsigned char a[n], b[m];
			FILE* file;
			errno_t err;
			err = fopen_s(&file, kontpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			FILE* msg;
			err = fopen_s(&msg, msgpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			remove("lsb.wav");
			FILE* filenew;
			err = fopen_s(&filenew, "lsb.wav", "wb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 0; i < kontdata/n; i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, kontdata%n, 1, file);
			fwrite(&a, kontdata%n, 1, filenew);
			for (int i = 1; i <= (msgsize / m); i++) {
				fread(&b, m, 1, msg);
				fread(&a, n, 1, file);
				long pt = 0;
				for (int j = 0; j < m; j++) {
					for (int k = 7; k >= 0; k--) {
						if (b[j] & (1 << k)) {
							a[pt] |= (1 << 0);
						}
						else {
							a[pt] &= ~(1 << 0);
						}
						pt += 2;
					}
				}
				fwrite(&a, n, 1, filenew);
			}
			fread(&b, msgsize % m, 1, msg);
			long ns;
			if ((long)(header.chunkSize - ((msgsize / m+ (((msgsize % m)>0)&& ((msgsize / m)==0))) * n+kontdata)) >= n) {
				ns = n;
			}
			else {
				ns = (header.chunkSize-kontdata) % n;
			}
			fread(&a, ns, 1, file);
			long pt = 0;
			for (int j = 0; j < msgsize % m; j++) {
				for (int k = 7; k >= 0; k--) {
					if (b[j] & (1 << k)) {
						a[pt] |= (1 << 0);
					}
					else {
						a[pt] &= ~(1 << 0);
					}
					pt += 2;
				}
			}
			fwrite(&a, ns, 1, filenew);
			if (ns == n) {
				
				for (int i = 1; i <= (long)((header.chunkSize-(kontdata+(msgsize/m+(msgsize%m>0))*n)) / n); i++) {
					fread(&a, n, 1, file);
					fwrite(&a, n, 1, filenew);
				}
				fread(&a, (long)((header.chunkSize - (kontdata + (msgsize / m + (msgsize%m > 0))*n)) % n), 1, file);
				fwrite(&a, (long)((header.chunkSize - (kontdata + (msgsize / m + (msgsize%m > 0))*n)) % n), 1, filenew);
			}
			fclose(filenew);
			fclose(msg);
			fclose(file);
		}

		void ParCode() 
		{
			unsigned char a[n], b[m];
			FILE* file;
			errno_t err;
			err = fopen_s(&file, kontpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			FILE* msg;
			err = fopen_s(&msg, msgpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			remove("paritycoding.wav");
			FILE* filenew;
			err = fopen_s(&filenew, "paritycoding.wav", "wb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 0; i < kontdata / n; i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, kontdata%n, 1, file);
			fwrite(&a, kontdata%n, 1, filenew);
			long blocksize = kontsize / (msgsize * 8);
			blocksize = (blocksize - blocksize % 2) * 8;
			long wr = kontdata;
			if (blocksize <= n) {
				long k = n / blocksize;
				for (int i = 0; i < msgsize / k; i++) {
					fread(&a, k*blocksize, 1, file);
					fread(&b, k, 1, msg);
					long pt = 0;
					for (int j = 0; j < k; j++) {
						for (int l = 7; l >= 0; l--) {
							long count = 0;
							for (int s = 0; s < blocksize/8; s++) {
								count += NumUnits(a[pt]);
								pt++;
							}
							bool bt = (bool((1 << k)&b[j]));
							if (count % 2 != bt) {
								a[pt - 2] ^= (1 << 0);
							}
						}
					}
					fwrite(&a, k*blocksize, 1, filenew);
				}
				fread(&a, (msgsize % k)*blocksize, 1, file);
				fread(&b, msgsize % k, 1, msg);
				long pt = 0;
				for (int j = 0; j < (msgsize % k); j++) {
					for (int l = 7; l >= 0; l--) {
						long count = 0;
						for (int s = 0; s < blocksize/8; s++) {
							count += NumUnits(a[pt]);
							pt++;
						}
						bool bt = (bool((1 << k)&b[j]));
						if (count % 2 != bt) {
							a[pt - 2] ^= (1 << 0);
						}
					}
				}
				fwrite(&a, (msgsize % k)*blocksize, 1, filenew);
				wr += blocksize * msgsize;
			}
			else {
				blocksize = blocksize / 8;
				if (blocksize <= n) {
					long z = blocksize;
					long t = 1;
					while (z * 2 <= n) {
						z = z * 2;
						t = t * 2;
					}
					for (int i = 0; i < msgsize; i++) {
						fread(&b, 1, 1, msg);
						for (int j = 0; j < 8 / t; j++) {
							fread(&a, z, 1, file);
							long pt = 0;
							for (int k = 0; k < t; k++) {
								long count = 0;
								for (int l = 0; l < blocksize; l++) {
									count += NumUnits(a[pt]);
									pt++;
								}
								bool bt = (bool((1 << (7-(j*t+k)))&b[0]));
								if (count % 2 != bt) {
									a[pt - 2] ^= (1 << 0);
								}
							}
							fwrite(&a, z, 1, filenew);
						}
					}
					wr+=blocksize*8*msgsize;
				}
				else {
					for (int i = 0; i < msgsize; i++) {
						fread(&b, 1, 1, msg);
						for (int j = 7; j >= 0; j--) {
							long long count = 0;
							for (int k = 0; k < blocksize / n; k++) {
								fread(&a, n, 1, file);
								for (int l = 0; l < n; l++) {
									count += NumUnits(a[l]);
								}
								fwrite(&a, n, 1, filenew);
							}
							fread(&a, blocksize % n, 1, file);
							for (int l = 0; l < blocksize % n; l++) {
								count += NumUnits(a[l]);
							}
							bool bt = (bool((1 << (j))&b[0]));
							if (count % 2 != bt) {
								a[blocksize % n - 2] ^= (1 << 0);
							}
							fwrite(&a, blocksize % n, 1, filenew);
						}
					}
					wr += blocksize * 8 * msgsize;
				}
			}
			for (int i = 0; i < (header.chunkSize - wr) / n; i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, (header.chunkSize - wr) % n, 1, file);
			fwrite(&a, (header.chunkSize - wr) % n, 1, filenew);
			fclose(file);
			fclose(filenew);
			fclose(msg);
		}

		void Echo()
		{
			short a[n/2+1000], c[500];
			unsigned char b[m-1000];
			FILE* file;
			errno_t err;
			err = fopen_s(&file, kontpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			FILE* msg;
			err = fopen_s(&msg, msgpath, "rb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			remove("echo.wav");
			FILE* filenew;
			err = fopen_s(&filenew, "echo.wav", "wb");
			if (err) {
				MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
				return;
			}
			for (int i = 0; i < kontdata / n; i++) {
				fread(&a, n, 1, file);
				fwrite(&a, n, 1, filenew);
			}
			fread(&a, kontdata%n, 1, file);
			fwrite(&a, kontdata%n, 1, filenew);
			bool last = 1;
			long k1 = long(0.001*byteRate);
			long k2 = long(0.0005*byteRate);
			long rk = n - (n % (byteRate / 2))+k1*2;
			long rm = rk / (byteRate / 2);
			const float base = 0.1;
			float kf = 0;
			long pt = k1;
			
			fread(&c, k1*2, 1, file);
			for (int i = 0; i < msgsize / rm; i++) {
				fread(&b, rm, 1, msg);
				fseek(file, -k1*2, SEEK_CUR);
				pt = k1;
				fread(&a, rk, 1, file);
				for (int j = 0; j < k1; j++) {
					a[j] = c[j];
				}
				for (int j = 0; j < rm; j++) {
					for (int x = 7; x >= 0; x--) {
						if (last == 1) {
							kf = base;
						}
						else kf = 0;
						last = (bool((1 << x)&b[j]));
						for (int k = pt; k < (j * 8 + (8 - x))*byteRate / 32+k1; k++) {
							
							short prov = a[k] + (short)(kf*a[k - k1]) + (short)((base - kf)*a[k - k2]);
							if (prov<32767 || prov>-32768) a[k] = prov;
							else {
								if (prov < 0) a[k] = -32768;
								else a[k] = 32767;
							}
							if (last) {
								if (kf < base) kf += base / (byteRate / 320);
								if (kf > base) kf = base;
							}
							else {
								if (kf > 0) kf -= base / (byteRate / 320);
								if (kf < 0) kf = 0;
							}
							pt++;
						}
					}
				}
				for (int k = 0; k < k1; k++) {
					c[k] = a[rk / 2 - k1 + k];
				}
				fwrite(&a, rk - k1 * 2, 1, filenew);
			}
			if (msgsize%rm > 0) {
				fread(&b, msgsize%rm, 1, msg);
				fseek(file, -k1 * 2, SEEK_CUR);
				pt = k1;
				fread(&a, (byteRate / 2)*(msgsize%rm)+k1*2, 1, file);
				for (int j = 0; j < k1; j++) {
					a[j] = c[j];
				}
				for (int j = 0; j < rm; j++) {
					for (int x = 7; x >= 0; x--) {
						if (last == 1) {
							kf = base;
						}
						else kf = 0;
						last = (bool((1 << x)&b[j]));
						for (int k = pt; k < (j * 8 + (8 - x))*byteRate / 32+k1; k++) {

							short prov = a[k] + (short)(kf*a[k - k1]) + (short)((base - kf)*a[k - k2]);
							if (prov<32767 || prov>-32768) a[k] = prov;
							else {
								if (prov < 0) a[k] = -32768;
								else a[k] = 32767;
							}
							if (last) {
								if (kf < base) kf += base / (byteRate / 320);
								if (kf > base) kf = base;
							}
							else {
								if (kf > 0) kf -= base / (byteRate / 320);
								if (kf < 0) kf = 0;
							}
							pt++;
						}
					}
				}
				for (int k = 0; k < k1; k++) {
					c[k] = a[((byteRate / 2)*(msgsize%rm))/2  + k];
				}
				fwrite(&a, (byteRate / 2)*(msgsize%rm), 1, filenew);
			}
			long rd = kontsize - (msgsize * byteRate / 2);
			pt = k1;
			if (last) { 
				kf = base; 
			}
			else kf = 0;
			for (int i = 0; i < rd / n; i++) {
				fseek(file, -k1 * 2, SEEK_CUR);
				fread(&a, n + k1 * 2, 1, file);
				for (int j = 0; j < k1; j++) {
					a[j] = c[j];
				}
				for (int j = pt; j < n/2 + k1; j++) {
					short prov = a[j] + (short)(kf*a[j - k1]) + (short)((base - kf)*a[j - k2]);
					if (prov<32767 || prov>-32768) a[j] = prov;
					else {
						if (prov < 0) a[j] = -32768;
						else a[j] = 32767;
					}
				}
				for (int k = 0; k < k1; k++) {
					c[k] = a[n/2 + k];
				}
				fwrite(&a, n , 1, filenew);
			}

			if (rd%n > 0) {
				fseek(file, -k1 * 2, SEEK_CUR);
				fread(&a, rd%n + k1 * 2, 1, file);
				for (int j = 0; j < k1; j++) {
					a[j] = c[j];
				}
				for (int j = pt; j < (rd%n) / 2 + k1; j++) {
					short prov = a[j] + (short)(kf*a[j - k1]) + (short)((base - kf)*a[j - k2]);
					if (prov<32767 || prov>-32768) a[j] = prov;
					else {
						if (prov < 0) a[j] = -32768;
						else a[j] = 32767;
					}
				}
				for (int k = 0; k < k1; k++) {
					c[k] = a[(rd%n) / 2 - k1 + k];
				}
				fwrite(&a, rd%n - k1 * 2, 1, filenew);
			}
			fwrite(&c, k1 * 2, 1, filenew);
			if (header.chunkSize - kontdata - kontsize > 0) {
				fread(&a, header.chunkSize - kontdata - kontsize, 1, file);
				fwrite(&a, header.chunkSize - kontdata - kontsize, 1, filenew);
			}
			fclose(file);
			fclose(filenew);
			fclose(msg);
		}
	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::OpenFileDialog^  openFileDialog1;
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Button^  button2;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::PictureBox^  pictureBox1;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::PictureBox^  pictureBox2;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::PictureBox^  pictureBox3;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::PictureBox^  pictureBox4;
private: System::Windows::Forms::Label^  label8;
private: System::Windows::Forms::PictureBox^  pictureBox5;
private: System::Windows::Forms::Label^  label9;
private: System::Windows::Forms::Label^  label10;
private: System::Windows::Forms::Label^  label11;
private: System::Windows::Forms::Label^  label12;
private: System::Windows::Forms::Label^  label13;
private: System::Windows::Forms::Label^  label14;
private: System::Windows::Forms::Label^  label15;
private: System::Windows::Forms::Label^  label16;
private: System::Windows::Forms::Label^  label17;
private: System::Windows::Forms::Label^  label18;
private: System::Windows::Forms::Label^  label19;
private: System::Windows::Forms::Label^  label20;
private: System::Windows::Forms::Label^  label21;
private: System::Windows::Forms::Label^  label22;
private: System::Windows::Forms::Label^  label23;
private: System::Windows::Forms::Label^  label24;
private: System::Windows::Forms::Label^  label25;
private: System::Windows::Forms::Label^  label26;
private: System::Windows::Forms::Label^  label27;
private: System::Windows::Forms::Label^  label28;
private: System::Windows::Forms::Label^  label29;
private: System::Windows::Forms::Label^  label30;
private: System::Windows::Forms::Label^  label31;
private: System::Windows::Forms::Label^  label32;
private: System::Windows::Forms::Label^  label33;
private: System::Windows::Forms::Label^  label34;
private: System::Windows::Forms::Label^  label35;
private: System::Windows::Forms::Label^  label36;


	protected:

	private:
		/// <summary>
		/// Обязательная переменная конструктора.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Требуемый метод для поддержки конструктора — не изменяйте 
		/// содержимое этого метода с помощью редактора кода.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(MyForm::typeid));
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->openFileDialog1 = (gcnew System::Windows::Forms::OpenFileDialog());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->pictureBox2 = (gcnew System::Windows::Forms::PictureBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->pictureBox3 = (gcnew System::Windows::Forms::PictureBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->pictureBox4 = (gcnew System::Windows::Forms::PictureBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->pictureBox5 = (gcnew System::Windows::Forms::PictureBox());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->label15 = (gcnew System::Windows::Forms::Label());
			this->label16 = (gcnew System::Windows::Forms::Label());
			this->label17 = (gcnew System::Windows::Forms::Label());
			this->label18 = (gcnew System::Windows::Forms::Label());
			this->label19 = (gcnew System::Windows::Forms::Label());
			this->label20 = (gcnew System::Windows::Forms::Label());
			this->label21 = (gcnew System::Windows::Forms::Label());
			this->label22 = (gcnew System::Windows::Forms::Label());
			this->label23 = (gcnew System::Windows::Forms::Label());
			this->label24 = (gcnew System::Windows::Forms::Label());
			this->label25 = (gcnew System::Windows::Forms::Label());
			this->label26 = (gcnew System::Windows::Forms::Label());
			this->label27 = (gcnew System::Windows::Forms::Label());
			this->label28 = (gcnew System::Windows::Forms::Label());
			this->label29 = (gcnew System::Windows::Forms::Label());
			this->label30 = (gcnew System::Windows::Forms::Label());
			this->label31 = (gcnew System::Windows::Forms::Label());
			this->label32 = (gcnew System::Windows::Forms::Label());
			this->label33 = (gcnew System::Windows::Forms::Label());
			this->label34 = (gcnew System::Windows::Forms::Label());
			this->label35 = (gcnew System::Windows::Forms::Label());
			this->label36 = (gcnew System::Windows::Forms::Label());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox4))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox5))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label1->Location = System::Drawing::Point(36, 25);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(199, 19);
			this->label1->TabIndex = 0;
			this->label1->Text = L"Выберите файл-контейнер:";
			this->label1->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// openFileDialog1
			// 
			this->openFileDialog1->FileName = L"Файл не выбран";
			// 
			// textBox1
			// 
			this->textBox1->Enabled = false;
			this->textBox1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox1->Location = System::Drawing::Point(233, 131);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(278, 27);
			this->textBox1->TabIndex = 1;
			this->textBox1->Text = L"Файл не выбран";
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button1->Location = System::Drawing::Point(40, 122);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(148, 42);
			this->button1->TabIndex = 2;
			this->button1->Text = L"Выбрать файл";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// button2
			// 
			this->button2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button2->Location = System::Drawing::Point(777, 122);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(148, 42);
			this->button2->TabIndex = 4;
			this->button2->Text = L"Выбрать файл";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// textBox2
			// 
			this->textBox2->Enabled = false;
			this->textBox2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->textBox2->Location = System::Drawing::Point(969, 131);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(278, 27);
			this->textBox2->TabIndex = 3;
			this->textBox2->Text = L"Файл не выбран";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label2->Location = System::Drawing::Point(773, 25);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(208, 19);
			this->label2->TabIndex = 5;
			this->label2->Text = L"Выберите файл-сообщение:";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label3->Location = System::Drawing::Point(773, 56);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(0, 19);
			this->label3->TabIndex = 6;
			// 
			// button3
			// 
			this->button3->Enabled = false;
			this->button3->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button3->Location = System::Drawing::Point(516, 207);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(312, 47);
			this->button3->TabIndex = 7;
			this->button3->Text = L"Начать процесс сокрытия данных";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &MyForm::button3_Click);
			// 
			// pictureBox1
			// 
			this->pictureBox1->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox1.Image")));
			this->pictureBox1->Location = System::Drawing::Point(40, 328);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(47, 50);
			this->pictureBox1->TabIndex = 8;
			this->pictureBox1->TabStop = false;
			this->pictureBox1->Visible = false;
			this->pictureBox1->Click += gcnew System::EventHandler(this, &MyForm::pictureBox1_Click);
			this->pictureBox1->DoubleClick += gcnew System::EventHandler(this, &MyForm::pictureBox1_Click);
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label4->Location = System::Drawing::Point(36, 286);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(145, 19);
			this->label4->TabIndex = 9;
			this->label4->Text = L"Загруженный файл:";
			this->label4->Visible = false;
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label5->Location = System::Drawing::Point(36, 397);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(123, 19);
			this->label5->TabIndex = 11;
			this->label5->Text = L"Метод оверлея:";
			this->label5->Visible = false;
			// 
			// pictureBox2
			// 
			this->pictureBox2->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox2.Image")));
			this->pictureBox2->Location = System::Drawing::Point(40, 441);
			this->pictureBox2->Name = L"pictureBox2";
			this->pictureBox2->Size = System::Drawing::Size(47, 50);
			this->pictureBox2->TabIndex = 10;
			this->pictureBox2->TabStop = false;
			this->pictureBox2->Visible = false;
			this->pictureBox2->Click += gcnew System::EventHandler(this, &MyForm::pictureBox2_Click);
			this->pictureBox2->DoubleClick += gcnew System::EventHandler(this, &MyForm::pictureBox2_Click);
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label6->Location = System::Drawing::Point(36, 513);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(206, 19);
			this->label6->TabIndex = 13;
			this->label6->Text = L"Метод Least Significant Bit:";
			this->label6->Visible = false;
			// 
			// pictureBox3
			// 
			this->pictureBox3->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox3.Image")));
			this->pictureBox3->Location = System::Drawing::Point(40, 560);
			this->pictureBox3->Name = L"pictureBox3";
			this->pictureBox3->Size = System::Drawing::Size(47, 50);
			this->pictureBox3->TabIndex = 12;
			this->pictureBox3->TabStop = false;
			this->pictureBox3->Visible = false;
			this->pictureBox3->Click += gcnew System::EventHandler(this, &MyForm::pictureBox3_Click);
			this->pictureBox3->DoubleClick += gcnew System::EventHandler(this, &MyForm::pictureBox3_Click);
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label7->Location = System::Drawing::Point(36, 634);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(216, 19);
			this->label7->TabIndex = 15;
			this->label7->Text = L"Метод четного кодирования:";
			this->label7->Visible = false;
			// 
			// pictureBox4
			// 
			this->pictureBox4->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox4.Image")));
			this->pictureBox4->Location = System::Drawing::Point(40, 680);
			this->pictureBox4->Name = L"pictureBox4";
			this->pictureBox4->Size = System::Drawing::Size(47, 50);
			this->pictureBox4->TabIndex = 14;
			this->pictureBox4->TabStop = false;
			this->pictureBox4->Visible = false;
			this->pictureBox4->Click += gcnew System::EventHandler(this, &MyForm::pictureBox4_Click);
			this->pictureBox4->DoubleClick += gcnew System::EventHandler(this, &MyForm::pictureBox4_Click);
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label8->Location = System::Drawing::Point(36, 754);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(91, 19);
			this->label8->TabIndex = 17;
			this->label8->Text = L"Эхо-метод:";
			this->label8->Visible = false;
			// 
			// pictureBox5
			// 
			this->pictureBox5->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox5.Image")));
			this->pictureBox5->Location = System::Drawing::Point(40, 800);
			this->pictureBox5->Name = L"pictureBox5";
			this->pictureBox5->Size = System::Drawing::Size(47, 50);
			this->pictureBox5->TabIndex = 16;
			this->pictureBox5->TabStop = false;
			this->pictureBox5->Visible = false;
			this->pictureBox5->Click += gcnew System::EventHandler(this, &MyForm::pictureBox5_Click);
			this->pictureBox5->DoubleClick += gcnew System::EventHandler(this, &MyForm::pictureBox5_Click);
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label9->Location = System::Drawing::Point(348, 286);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(110, 19);
			this->label9->TabIndex = 18;
			this->label9->Text = L"Размер файла:";
			this->label9->Visible = false;
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label10->Location = System::Drawing::Point(348, 397);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(197, 19);
			this->label10->TabIndex = 19;
			this->label10->Text = L"Изменение размера файла:";
			this->label10->Visible = false;
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label11->Location = System::Drawing::Point(348, 513);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(197, 19);
			this->label11->TabIndex = 20;
			this->label11->Text = L"Изменение размера файла:";
			this->label11->Visible = false;
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label12->Location = System::Drawing::Point(348, 634);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(197, 19);
			this->label12->TabIndex = 21;
			this->label12->Text = L"Изменение размера файла:";
			this->label12->Visible = false;
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label13->Location = System::Drawing::Point(348, 754);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(197, 19);
			this->label13->TabIndex = 22;
			this->label13->Text = L"Изменение размера файла:";
			this->label13->Visible = false;
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label14->Location = System::Drawing::Point(348, 328);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(0, 19);
			this->label14->TabIndex = 23;
			this->label14->Visible = false;
			// 
			// label15
			// 
			this->label15->AutoSize = true;
			this->label15->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label15->Location = System::Drawing::Point(348, 441);
			this->label15->Name = L"label15";
			this->label15->Size = System::Drawing::Size(0, 19);
			this->label15->TabIndex = 24;
			this->label15->Visible = false;
			// 
			// label16
			// 
			this->label16->AutoSize = true;
			this->label16->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label16->Location = System::Drawing::Point(348, 560);
			this->label16->Name = L"label16";
			this->label16->Size = System::Drawing::Size(0, 19);
			this->label16->TabIndex = 25;
			this->label16->Visible = false;
			// 
			// label17
			// 
			this->label17->AutoSize = true;
			this->label17->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label17->Location = System::Drawing::Point(348, 680);
			this->label17->Name = L"label17";
			this->label17->Size = System::Drawing::Size(0, 19);
			this->label17->TabIndex = 26;
			this->label17->Visible = false;
			// 
			// label18
			// 
			this->label18->AutoSize = true;
			this->label18->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label18->Location = System::Drawing::Point(348, 800);
			this->label18->Name = L"label18";
			this->label18->Size = System::Drawing::Size(0, 19);
			this->label18->TabIndex = 27;
			this->label18->Visible = false;
			// 
			// label19
			// 
			this->label19->AutoSize = true;
			this->label19->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label19->Location = System::Drawing::Point(685, 397);
			this->label19->Name = L"label19";
			this->label19->Size = System::Drawing::Size(247, 19);
			this->label19->TabIndex = 28;
			this->label19->Text = L"Среднеквадратичное отклонение:";
			this->label19->Visible = false;
			// 
			// label20
			// 
			this->label20->AutoSize = true;
			this->label20->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label20->Location = System::Drawing::Point(685, 513);
			this->label20->Name = L"label20";
			this->label20->Size = System::Drawing::Size(247, 19);
			this->label20->TabIndex = 29;
			this->label20->Text = L"Среднеквадратичное отклонение:";
			this->label20->Visible = false;
			// 
			// label21
			// 
			this->label21->AutoSize = true;
			this->label21->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label21->Location = System::Drawing::Point(685, 634);
			this->label21->Name = L"label21";
			this->label21->Size = System::Drawing::Size(247, 19);
			this->label21->TabIndex = 30;
			this->label21->Text = L"Среднеквадратичное отклонение:";
			this->label21->Visible = false;
			// 
			// label22
			// 
			this->label22->AutoSize = true;
			this->label22->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label22->Location = System::Drawing::Point(685, 754);
			this->label22->Name = L"label22";
			this->label22->Size = System::Drawing::Size(247, 19);
			this->label22->TabIndex = 31;
			this->label22->Text = L"Среднеквадратичное отклонение:";
			this->label22->Visible = false;
			// 
			// label23
			// 
			this->label23->AutoSize = true;
			this->label23->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label23->Location = System::Drawing::Point(685, 441);
			this->label23->Name = L"label23";
			this->label23->Size = System::Drawing::Size(0, 19);
			this->label23->TabIndex = 32;
			this->label23->Visible = false;
			// 
			// label24
			// 
			this->label24->AutoSize = true;
			this->label24->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label24->Location = System::Drawing::Point(685, 560);
			this->label24->Name = L"label24";
			this->label24->Size = System::Drawing::Size(0, 19);
			this->label24->TabIndex = 33;
			this->label24->Visible = false;
			// 
			// label25
			// 
			this->label25->AutoSize = true;
			this->label25->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label25->Location = System::Drawing::Point(685, 680);
			this->label25->Name = L"label25";
			this->label25->Size = System::Drawing::Size(0, 19);
			this->label25->TabIndex = 34;
			this->label25->Visible = false;
			// 
			// label26
			// 
			this->label26->AutoSize = true;
			this->label26->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label26->Location = System::Drawing::Point(685, 800);
			this->label26->Name = L"label26";
			this->label26->Size = System::Drawing::Size(0, 19);
			this->label26->TabIndex = 35;
			this->label26->Visible = false;
			// 
			// label27
			// 
			this->label27->AutoSize = true;
			this->label27->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label27->Location = System::Drawing::Point(773, 89);
			this->label27->Name = L"label27";
			this->label27->Size = System::Drawing::Size(0, 19);
			this->label27->TabIndex = 36;
			// 
			// label28
			// 
			this->label28->AutoSize = true;
			this->label28->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label28->Location = System::Drawing::Point(36, 56);
			this->label28->Name = L"label28";
			this->label28->Size = System::Drawing::Size(416, 19);
			this->label28->TabIndex = 37;
			this->label28->Text = L"(для применения эхо-метода необходим файл без сжатия)";
			// 
			// label29
			// 
			this->label29->AutoSize = true;
			this->label29->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label29->Location = System::Drawing::Point(1077, 800);
			this->label29->Name = L"label29";
			this->label29->Size = System::Drawing::Size(0, 19);
			this->label29->TabIndex = 45;
			this->label29->Visible = false;
			// 
			// label30
			// 
			this->label30->AutoSize = true;
			this->label30->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label30->Location = System::Drawing::Point(1077, 680);
			this->label30->Name = L"label30";
			this->label30->Size = System::Drawing::Size(0, 19);
			this->label30->TabIndex = 44;
			this->label30->Visible = false;
			// 
			// label31
			// 
			this->label31->AutoSize = true;
			this->label31->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label31->Location = System::Drawing::Point(1077, 560);
			this->label31->Name = L"label31";
			this->label31->Size = System::Drawing::Size(0, 19);
			this->label31->TabIndex = 43;
			this->label31->Visible = false;
			// 
			// label32
			// 
			this->label32->AutoSize = true;
			this->label32->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label32->Location = System::Drawing::Point(1077, 441);
			this->label32->Name = L"label32";
			this->label32->Size = System::Drawing::Size(0, 19);
			this->label32->TabIndex = 42;
			this->label32->Visible = false;
			// 
			// label33
			// 
			this->label33->AutoSize = true;
			this->label33->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label33->Location = System::Drawing::Point(1077, 754);
			this->label33->Name = L"label33";
			this->label33->Size = System::Drawing::Size(270, 19);
			this->label33->TabIndex = 41;
			this->label33->Text = L"Пиковое отношение сигнала к шуму:";
			this->label33->Visible = false;
			// 
			// label34
			// 
			this->label34->AutoSize = true;
			this->label34->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label34->Location = System::Drawing::Point(1077, 634);
			this->label34->Name = L"label34";
			this->label34->Size = System::Drawing::Size(270, 19);
			this->label34->TabIndex = 40;
			this->label34->Text = L"Пиковое отношение сигнала к шуму:";
			this->label34->Visible = false;
			// 
			// label35
			// 
			this->label35->AutoSize = true;
			this->label35->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label35->Location = System::Drawing::Point(1077, 513);
			this->label35->Name = L"label35";
			this->label35->Size = System::Drawing::Size(270, 19);
			this->label35->TabIndex = 39;
			this->label35->Text = L"Пиковое отношение сигнала к шуму:";
			this->label35->Visible = false;
			// 
			// label36
			// 
			this->label36->AutoSize = true;
			this->label36->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label36->Location = System::Drawing::Point(1077, 397);
			this->label36->Name = L"label36";
			this->label36->Size = System::Drawing::Size(270, 19);
			this->label36->TabIndex = 38;
			this->label36->Text = L"Пиковое отношение сигнала к шуму:";
			this->label36->Visible = false;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackColor = System::Drawing::SystemColors::ControlLightLight;
			this->ClientSize = System::Drawing::Size(1384, 862);
			this->Controls->Add(this->label29);
			this->Controls->Add(this->label30);
			this->Controls->Add(this->label31);
			this->Controls->Add(this->label32);
			this->Controls->Add(this->label33);
			this->Controls->Add(this->label34);
			this->Controls->Add(this->label35);
			this->Controls->Add(this->label36);
			this->Controls->Add(this->label28);
			this->Controls->Add(this->label27);
			this->Controls->Add(this->label26);
			this->Controls->Add(this->label25);
			this->Controls->Add(this->label24);
			this->Controls->Add(this->label23);
			this->Controls->Add(this->label22);
			this->Controls->Add(this->label21);
			this->Controls->Add(this->label20);
			this->Controls->Add(this->label19);
			this->Controls->Add(this->label18);
			this->Controls->Add(this->label17);
			this->Controls->Add(this->label16);
			this->Controls->Add(this->label15);
			this->Controls->Add(this->label14);
			this->Controls->Add(this->label13);
			this->Controls->Add(this->label12);
			this->Controls->Add(this->label11);
			this->Controls->Add(this->label10);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->pictureBox5);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->pictureBox4);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->pictureBox3);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->pictureBox2);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->FormBorderStyle = System::Windows::Forms::FormBorderStyle::FixedSingle;
			this->MaximizeBox = false;
			this->Name = L"MyForm";
			this->StartPosition = System::Windows::Forms::FormStartPosition::CenterScreen;
			this->Text = L"Изучение методов стеганографического сокрытия информации";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load_1);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox3))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox4))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox5))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void label1_Click(System::Object^  sender, System::EventArgs^  e) {
		
	}
			
private: System::Void pictureBox1_Click(System::Object^  sender, System::EventArgs^  e) {
		PlaySound(NULL, 0, 0);
		bool l = nowplaying[0];
		for (int i = 0; i < 5; i++) {
			nowplaying[i] = false;
		}
		pictureBox1->Image = Image::FromFile("play.png");
		pictureBox2->Image = Image::FromFile("play.png");
		pictureBox3->Image = Image::FromFile("play.png");
		pictureBox4->Image = Image::FromFile("play.png");
		pictureBox5->Image = Image::FromFile("play.png");
		if (l == false) {
			nowplaying[0] = true;
			PlaySound(L"orig.wav", NULL, SND_ASYNC| SND_FILENAME | SND_LOOP);
			pictureBox1->Image = Image::FromFile("stop.png");
		}
	}
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
	openFileDialog1->Filter = "All files (*.*)|*.*";
	if (openFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK)
	{
		msgchecked = false;
		button3->Enabled = (kontchecked&&msgchecked);
		msgsize = 0;
		long len = openFileDialog1->FileName->Length;
		for (int i = 0; i < len; i++) {
			msgpath[i] = openFileDialog1->FileName[i];
		}
		for (int i = len; i < 256; i++) {
			msgpath[i] = NULL;
		}
		FILE* file;
		errno_t err;
		err = fopen_s(&file, msgpath, "rb");
		if (err) {
			MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
			return;
		}
		fclose(file);
		fstream msg(msgpath);
		msg.seekg(0, ios::end);
		msgsize = msg.tellg();
		msg.close();
		textBox2->Text = openFileDialog1->FileName;
		msgchecked = true;
		button3->Enabled = (kontchecked&&msgchecked);

	}
}
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
	openFileDialog1->Filter = "WAVE files (*.wav)|*.wav";
	if (openFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK)
	{
		textBox1->Text = "Файл не выбран";
		label3->Text = "";
		label27->Text = "";
		kontchecked = false;
		button3->Enabled = (kontchecked&&msgchecked);
		kontdata = 0;
		kontsize = 0;
		header.chunkId = 0;
		header.chunkSize = 0;
		header.format = 0;
		long len = openFileDialog1->FileName->Length;
		for (int i = 0; i < len; i++) {
			kontpath[i] = openFileDialog1->FileName[i];
		}
		for (int i = len; i < 256; i++) {
			kontpath[i] = NULL;
		}
		FILE* file;
		errno_t err;
		err = fopen_s(&file, kontpath, "rb");
		if (err) {
			MessageBox::Show("Возникла ошибка при открытии файла!", "Ошибка");
			return;
		}
		fseek(file, 20, SEEK_SET);
		fread_s(&audioFormat, sizeof(short), sizeof(short), 1, file);
		fseek(file, 28, SEEK_SET);
		fread_s(&byteRate, sizeof(long), sizeof(long), 1, file);
		fseek(file, 0, SEEK_SET);
		fread_s(&header, sizeof(WAVHEADER), sizeof(WAVHEADER), 1, file);
		if (header.chunkId != 0x46464952 || header.format != 0x45564157) {
			MessageBox::Show("Возникла ошибка при чтении файла!", "Ошибка");
			return;
		}
		header.chunkSize += 8;
		long long pt = 12;
		long chunk = 0, chsize = 0;
		fread(&chunk, sizeof(long), 1, file);
		fread(&chsize, sizeof(long), 1, file);
		while (chunk != 0x61746164 && pt <= header.chunkSize) {
			pt += chsize + 8;
			fseek(file, pt, SEEK_SET);
			fread(&chunk, sizeof(long), 1, file);
			fread(&chsize, sizeof(long), 1, file);
		}
		if (chunk != 0x61746164) {
			MessageBox::Show("Возникла ошибка при чтении файла!", "Ошибка");
			return;
		}
		kontdata = pt + 8;
		kontsize = chsize;
		textBox1->Text = openFileDialog1->FileName;
		if (kontsize / 16>50000){
			label3->Text = "Максимальный размер сообщения: " + Convert::ToString(kontsize / 16/1024) + " Кб";
		}
		else label3->Text = "Максимальный размер сообщения: " + Convert::ToString(kontsize / 16) + " байт";
		if (audioFormat==1)
		label27->Text = " (для применения эхо-метода максимальный размер сообщения - " + Convert::ToString(long(kontsize / (byteRate / 2))) + " байт)";
		else MessageBox::Show("Данный файл имеет сжатие, эхо-метод применен не будет.", "Предупреждение");
		kontchecked = true;
		button3->Enabled = (kontchecked&&msgchecked);
	}
}
private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) {
	PlaySound(NULL, 0, 0);
	pictureBox1->Image = Image::FromFile("play.png");
	pictureBox2->Image = Image::FromFile("play.png");
	pictureBox3->Image = Image::FromFile("play.png");
	pictureBox4->Image = Image::FromFile("play.png");
	pictureBox5->Image = Image::FromFile("play.png");
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	if (msgsize > kontsize / 16) {
		MessageBox::Show("Сообщение превышает максимально допустимый объем!", "Ошибка");
		return;
	}
	label4->Visible = false;
	pictureBox1->Visible = false;
	label5->Visible = false;
	pictureBox2->Visible = false;
	label6->Visible = false;
	pictureBox3->Visible = false;
	label7->Visible = false;
	pictureBox4->Visible = false;
	label8->Visible = false;
	pictureBox5->Visible = false;
	label9->Visible = false;
	label10->Visible = false;
	label11->Visible = false;
	label12->Visible = false;
	label13->Visible = false;
	label14->Text = Convert::ToString(long long(header.chunkSize) / 1024) + " Кб";
	label15->Text = Convert::ToString(msgsize) + " байт";
	label16->Text = "Нет";
	label17->Text = "Нет";
	label18->Text = "Нет";
	label18->Text = "Нет";
	char orig[256] = "orig.wav";
	for (int i = 8; i < 256; i++) {
		orig[i] = NULL;
	}
	label14->Visible = false;
	label15->Visible = false;
	label16->Visible = false;
	label17->Visible = false;
	label18->Visible = false;
	label19->Visible = false;
	label20->Visible = false;
	label21->Visible = false;
	label22->Visible = false;
	label23->Visible = false;
	label24->Visible = false;
	label25->Visible = false;
	label26->Visible = false;
	label29->Visible = false;
	label30->Visible = false;
	label31->Visible = false;
	label32->Visible = false;
	label33->Visible = false;
	label34->Visible = false;
	label35->Visible = false;
	label36->Visible = false;
	FileCopy();
	label4->Visible = true;
	pictureBox1->Visible = true;
	label9->Visible = true;
	label14->Visible = true;
	AddEnd();
	double psnr;
	double nrms = StandDev(orig, "addend.wav", psnr);
	label23->Text = Convert::ToString(nrms);
	label32->Text = Convert::ToString("-");
	label5->Visible = true;
	pictureBox2->Visible = true;
	label10->Visible = true;
	label15->Visible = true;
	label19->Visible = true;
	label23->Visible = true;
	label36->Visible = true;
	label32->Visible = true;
	LSB();
	nrms = StandDev(orig, "lsb.wav", psnr);
	psnr = long(20 * log10(psnr / nrms));
	label24->Text = Convert::ToString(nrms);
	label31->Text = Convert::ToString(psnr + " Дб");
	label6->Visible = true;
	pictureBox3->Visible = true;
	label11->Visible = true;
	label16->Visible = true;
	label20->Visible = true;
	label24->Visible = true;
	label35->Visible = true;
	label31->Visible = true;
	ParCode();
	nrms = StandDev(orig, "paritycoding.wav", psnr);
	psnr = long(20 * log10(psnr / nrms));
	label25->Text = Convert::ToString(nrms);
	label30->Text = Convert::ToString(psnr + " Дб");
	label7->Visible = true;
	pictureBox4->Visible = true;
	label12->Visible = true;
	label17->Visible = true;
	label21->Visible = true;
	label25->Visible = true;
	label34->Visible = true;
	label30->Visible = true;
	if (audioFormat == 1 && msgsize < kontsize / (byteRate / 2)) {
		Echo();
		nrms = StandDev(orig, "echo.wav", psnr);
		psnr = long(20 * log10(psnr / nrms));
		label26->Text = Convert::ToString(nrms);
		label29->Text = Convert::ToString(psnr + " Дб");
		label8->Visible = true;
		pictureBox5->Visible = true;
		label13->Visible = true;
		label18->Visible = true;
		label22->Visible = true;
		label26->Visible = true;
		label33->Visible = true;
		label29->Visible = true;
	}
	else remove("echo.wav");
}
private: System::Void pictureBox2_Click(System::Object^  sender, System::EventArgs^  e) {
	PlaySound(NULL, 0, 0);
	bool l = nowplaying[1];
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	pictureBox1->Image = Image::FromFile("play.png");
	pictureBox2->Image = Image::FromFile("play.png");
	pictureBox3->Image = Image::FromFile("play.png");
	pictureBox4->Image = Image::FromFile("play.png");
	pictureBox5->Image = Image::FromFile("play.png");
	if (l == false) {
		nowplaying[1] = true;
		PlaySound(L"addend.wav", NULL, SND_ASYNC | SND_FILENAME | SND_LOOP);
		pictureBox2->Image = Image::FromFile("stop.png");
	}
}
private: System::Void pictureBox3_Click(System::Object^  sender, System::EventArgs^  e) {
	PlaySound(NULL, 0, 0);
	bool l = nowplaying[2];
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	pictureBox1->Image = Image::FromFile("play.png");
	pictureBox2->Image = Image::FromFile("play.png");
	pictureBox3->Image = Image::FromFile("play.png");
	pictureBox4->Image = Image::FromFile("play.png");
	pictureBox5->Image = Image::FromFile("play.png");
	if (l == false) {
		nowplaying[2] = true;
		PlaySound(L"lsb.wav", NULL, SND_ASYNC | SND_FILENAME | SND_LOOP);
		pictureBox3->Image = Image::FromFile("stop.png");
	}
}
private: System::Void pictureBox4_Click(System::Object^  sender, System::EventArgs^  e) {
	PlaySound(NULL, 0, 0);
	bool l = nowplaying[3];
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	pictureBox1->Image = Image::FromFile("play.png");
	pictureBox2->Image = Image::FromFile("play.png");
	pictureBox3->Image = Image::FromFile("play.png");
	pictureBox4->Image = Image::FromFile("play.png");
	pictureBox5->Image = Image::FromFile("play.png");
	if (l == false) {
		nowplaying[3] = true;
		PlaySound(L"paritycoding.wav", NULL, SND_ASYNC | SND_FILENAME | SND_LOOP);
		pictureBox4->Image = Image::FromFile("stop.png");
	}
}
private: System::Void pictureBox5_Click(System::Object^  sender, System::EventArgs^  e) {
	PlaySound(NULL, 0, 0);
	bool l = nowplaying[4];
	for (int i = 0; i < 5; i++) {
		nowplaying[i] = false;
	}
	pictureBox1->Image = Image::FromFile("play.png");
	pictureBox2->Image = Image::FromFile("play.png");
	pictureBox3->Image = Image::FromFile("play.png");
	pictureBox4->Image = Image::FromFile("play.png");
	pictureBox5->Image = Image::FromFile("play.png");
	if (l == false) {
		nowplaying[4] = true;
		PlaySound(L"echo.wav", NULL, SND_ASYNC | SND_FILENAME | SND_LOOP);
		pictureBox5->Image = Image::FromFile("stop.png");
	}
}
private: System::Void MyForm_Load(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void MyForm_Load_1(System::Object^  sender, System::EventArgs^  e) {
}
};
}
