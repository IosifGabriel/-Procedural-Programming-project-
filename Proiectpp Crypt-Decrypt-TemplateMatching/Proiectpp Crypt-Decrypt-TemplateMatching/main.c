#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable : 4996)  // a se sterge daca nu e rulat pe Visual Studio


typedef struct
{
	double correlation;
	int i, j, ok;
	unsigned char culoare[3];


}fer;


void grayscale(unsigned char* MatLin)
{
	int Width = MatLin[19] * 256 + MatLin[18];
	int Heigth = MatLin[23] * 256 + MatLin[22];
	int i;
	for (i = 0; i < Width*Heigth; i++)
	{
		unsigned char aux = 0.299*MatLin[i * 3 + 56] + 0.587*MatLin[i * 3 + 55] + 0.114*MatLin[i * 3 + 54];
		MatLin[i * 3 + 56] = MatLin[i * 3 + 55] = MatLin[i * 3 + 54] = aux;

	}
}


unsigned int* XORSHIFT32(unsigned int seed, unsigned int lungime)
{
	unsigned int *R;
	R = (unsigned int*)malloc(sizeof(unsigned int)*(2 * lungime + 1));
	unsigned int k, r;
	R[0] = seed;
	r = seed;
	for (k = 1; k <= 2 * lungime; k++)
	{
		r = r ^ r << 13;
		r = r ^ r >> 17;
		r = r ^ r << 5;
		R[k] = r;
	}

	return R;
}


unsigned char* ReadBmpImg(char *ImgBmp)
{
	FILE * file;
	file = fopen(ImgBmp, "rb");

	if (file == NULL)
	{
		printf("eroare la deschiderea fisierului");
		return 0;
	}
	else
	{
		int valbit;
		int header[55];
		int index = 0;

		do
		{
			valbit = fgetc(file);
			header[index++] = valbit;
			if (index == 54)
				break;

		} while (valbit != EOF);

		int Width, Heigth;   // width-coloane heigth-linii
		Width = header[19] * 256 + header[18];
		Heigth = header[23] * 256 + header[22];
		int padding = (4 - ((Width * 3) % 4)) % 4;
		unsigned char **a;  // aloc si citesc matricea
		int i, j;

		a = (unsigned char**)malloc((Heigth) * sizeof(unsigned char*));
		for (i = 0; i < Heigth; i++)
		{
			*(a + i) = (unsigned char*)malloc((Width * 3 + padding) * sizeof(unsigned char));
		}

		fseek(file, 54, SEEK_SET);

		for (i = 0; i < Heigth; i++)
		{
			for (j = 0; j < Width * 3 + padding; j++)
				fread(&*(*(a + i) + j), sizeof(unsigned char), sizeof(unsigned char), file);
		}


		unsigned char *MatLin;   // am declarat matricea liniarizata si am citit-o
		MatLin = (unsigned char*)malloc(((Heigth * Width * 3) + 54) * sizeof(unsigned char));

		for (i = 0; i < 54; i++)
			MatLin[i] = header[i];


		for (i = Heigth - 1; i >= 0; i--)
			for (j = 0; j < Width * 3; j++)
				MatLin[(Heigth - 1 - i)*(Width * 3) + j + 54] = *(*(a + i) + j);

		for (i = 0; i < Heigth; i++)
			free(a[i]);
		free(a);


		return MatLin;
	}

	fclose(file);
}



void WriteBmpImg(char *BmpImg2, unsigned char *MatLinCrypt)
{
	FILE *fout;
	fout = fopen(BmpImg2, "wb");
	rewind(fout);
	int Width = MatLinCrypt[19] * 256 + MatLinCrypt[18];
	int Heigth = MatLinCrypt[23] * 256 + MatLinCrypt[22];
	int padding = (4 - ((Width * 3) % 4)) % 4;
	fwrite(MatLinCrypt, 1, 54, fout);
	int i;

	for (i = Heigth - 1; i >= 0; i--)
	{
		unsigned char pad[4] = { 0x00,0x00,0x00,0x00 };
		fwrite(MatLinCrypt + i * (Width * 3) + 54, sizeof(unsigned char), Width * 3, fout);
		if (padding != 0)
			fwrite(pad, sizeof(unsigned char), padding, fout);

	}

	fclose(fout);
}


void chipatrat(unsigned char *MatLinCrypt)
{
	int Width = MatLinCrypt[19] * 256 + MatLinCrypt[18];
	int Heigth = MatLinCrypt[23] * 256 + MatLinCrypt[22];

	unsigned double ficubara = Width * Heigth / 256;
	double xciudatrosu = 0;
	double xciudatverde = 0;
	double xciudatalbastru = 0;
	int i;


	for (i = 0; i < 256; i++)
	{
		int j = 0;
		double frecvrosu = 0, frecvverde = 0, frecvalbastru = 0;

		for (j = 56; j < Width*Heigth * 3; j += 3)
			if (i == MatLinCrypt[j])
				frecvrosu++;

		for (j = 55; j < Width*Heigth * 3; j += 3)
			if (i == MatLinCrypt[j])
				frecvverde++;

		for (j = 54; j < Width*Heigth * 3; j += 3)
			if (i == MatLinCrypt[j])
				frecvalbastru++;


		xciudatrosu += ((frecvrosu - ficubara)*(frecvrosu - ficubara)) / ficubara;
		xciudatverde += ((frecvverde - ficubara)*(frecvverde - ficubara)) / ficubara;
		xciudatalbastru += ((frecvalbastru - ficubara)*(frecvalbastru - ficubara)) / ficubara;

	}

	printf("CANAL ROSU: %.2lf \n", xciudatrosu);
	printf("CANAL VERDE: %.2lf \n", xciudatverde);
	printf("CANAL ALBASTRU: %.2lf \n", xciudatalbastru);

}


void Crypting(char *ImgBmp, char *ImgBmp2, char *SecretKey)
{
	unsigned char *MatLin = ReadBmpImg(ImgBmp);
	int Width = MatLin[19] * 256 + MatLin[18];
	int Heigth = MatLin[23] * 256 + MatLin[22];
	int padding = (4 - ((Width * 3) % 4)) % 4;
	FILE *fileseed = fopen(SecretKey, "r");
	unsigned  int seed, SV;
	fscanf(fileseed, "%u %u", &seed, &SV);
	unsigned  int *R = XORSHIFT32(seed, Width*Heigth);

	unsigned int *Perm;
	Perm = (unsigned int*)malloc(sizeof(unsigned int)*Width*Heigth);

	int i;
	for (i = 0; i < Width*Heigth; i++)
		Perm[i] = i;

	int j = 1;
	for (i = Heigth * Width - 1; i >= 1; i--)
	{
		unsigned  int indice = R[j++] % (i + 1);
		unsigned  int aux;
		aux = Perm[indice];
		Perm[indice] = Perm[i];
		Perm[i] = aux;

	}

	unsigned char *MatLinPerm;
	MatLinPerm = (unsigned char*)malloc(((Heigth * Width * 3) + 54) * sizeof(unsigned char));

	for (i = Width * Heigth - 1; i >= 0; i--)
	{
		MatLinPerm[Perm[i] * 3 + 54] = MatLin[i * 3 + 54];
		MatLinPerm[Perm[i] * 3 + 55] = MatLin[i * 3 + 55];
		MatLinPerm[Perm[i] * 3 + 56] = MatLin[i * 3 + 56];

	}

	j = Width * Heigth;

	unsigned char byteSV[3], byte[3];

	byteSV[2] = (SV >> 16) & 0xFF;
	byteSV[1] = (SV >> 8) & 0xFF;
	byteSV[0] = SV & 0xFF;

	byte[2] = (R[j] >> 16) & 0xFF;
	byte[1] = (R[j] >> 8) & 0xFF;
	byte[0] = R[j++] & 0xFF;

	unsigned char *MatLinCrypt;
	MatLinCrypt = (unsigned char*)malloc(sizeof(unsigned char)*(Width*Heigth * 3 + 54));

	for (i = 0; i < 54; i++)
		MatLinCrypt[i] = MatLin[i];

	
	MatLinCrypt[54] = byteSV[0] ^ MatLinPerm[54] ^ byte[0];
	MatLinCrypt[54 + 1] = byteSV[0 + 1] ^ MatLinPerm[54 + 1] ^ byte[0 + 1];
	MatLinCrypt[54 + 2] = byteSV[0 + 2] ^ MatLinPerm[54 + 2] ^ byte[0 + 2];
	
	for (i = 57; i < Width*Heigth * 3 + 54; i += 3)
	{
		byte[2] = (R[j] >> 16) & 0xFF;
		byte[1] = (R[j] >> 8) & 0xFF;
		byte[0] = R[j++] & 0xFF;

		MatLinCrypt[i] = MatLinCrypt[i - 3] ^ MatLinPerm[i] ^ byte[0];
		MatLinCrypt[i + 1] = MatLinCrypt[i - 3 + 1] ^ MatLinPerm[i + 1] ^ byte[0 + 1];
		MatLinCrypt[i + 2] = MatLinCrypt[i - 3 + 2] ^ MatLinPerm[i + 2] ^ byte[0 + 2];
	}
	printf("Test Chi-Squared imagine initiala: \n");
	chipatrat(MatLin);
	printf("\n");
	free(MatLin);
	free(R);
	free(MatLinPerm);
	WriteBmpImg(ImgBmp2, MatLinCrypt);
	printf("Test Chi-Squared imagine Criptata: \n");
	chipatrat(MatLinCrypt);
	free(MatLinCrypt);
	fclose(fileseed);

}


void Decrypt(char* BmpImg2, char* BmpImg3, char* SecretKey)
{
	unsigned char *MatLinCrypt = ReadBmpImg(BmpImg2);
	int Width = MatLinCrypt[19] * 256 + MatLinCrypt[18];
	int Heigth = MatLinCrypt[23] * 256 + MatLinCrypt[22];
	int padding = (4 - ((Width * 3) % 4)) % 4;

	FILE *fileseed = fopen(SecretKey, "r");
	unsigned  int seed, SV;
	fscanf(fileseed, "%u %u", &seed, &SV);
	unsigned  int *R;
	R = XORSHIFT32(seed, Width*Heigth);

	unsigned int *Perm;
	Perm = (unsigned int*)malloc(sizeof(unsigned int)*(Width * Heigth));

	int i;
	for (i = 0; i < Width*Heigth; i++)
		Perm[i] = i;

	int j = 1;
	for (i = Width * Heigth - 1; i >= 1; i--)
	{
		unsigned int indice = R[j++] % (i + 1);

		unsigned int aux;
		aux = Perm[i];
		Perm[i] = Perm[indice];
		Perm[indice] = aux;
	}

	unsigned  int *PermInv;
	PermInv = (unsigned  int*)malloc(sizeof(unsigned int)* (Width * Heigth));

	for (i = 0; i < Width*Heigth; i++)
		PermInv[Perm[i]]=i;

	unsigned char *MatLinCryptPrim;
	MatLinCryptPrim = (unsigned char*)malloc(sizeof(unsigned char)*((Width*Heigth * 3) + 54));

	unsigned char byteSV[3], byte[3];

	byteSV[2] = (SV >> 16) & 0xFF;
	byteSV[1] = (SV >> 8) & 0xFF;
	byteSV[0] = SV & 0xFF;
	j = Width * Heigth;
	byte[2] = (R[j] >> 16) & 0xFF;
	byte[1] = (R[j] >> 8) & 0xFF;
	byte[0] = R[j++] & 0xFF;


	MatLinCryptPrim[54] = byteSV[0] ^ MatLinCrypt[54] ^ byte[0];
	MatLinCryptPrim[54 + 1] = byteSV[0 + 1] ^ MatLinCrypt[54 + 1] ^ byte[0 + 1];
	MatLinCryptPrim[54 + 2] = byteSV[0 + 2] ^ MatLinCrypt[54 + 2] ^ byte[0 + 2];

	for (i = 57; i < Width*Heigth * 3 +54; i += 3)
	{
		byte[2] = (R[j] >> 16) & 0xFF;
		byte[1] = (R[j] >> 8) & 0xFF;
		byte[0] = R[j++] & 0xFF;

		MatLinCryptPrim[i] = MatLinCrypt[i - 3] ^ MatLinCrypt[i] ^ byte[0];
		MatLinCryptPrim[i + 1] = MatLinCrypt[i - 3 + 1] ^ MatLinCrypt[i + 1] ^ byte[0 + 1];
		MatLinCryptPrim[i + 2] = MatLinCrypt[i - 3 + 2] ^ MatLinCrypt[i + 2] ^ byte[0 + 2];
	}
	
	unsigned char *MatLinDecrypt;
	MatLinDecrypt = (unsigned char*)malloc(sizeof(unsigned char)*((Width*Heigth * 3) + 54));

	for (i = 0; i < Width*Heigth; i++)
	{
		MatLinDecrypt[PermInv[i] * 3 + 54] = MatLinCryptPrim[i * 3 + 54];
		MatLinDecrypt[PermInv[i] * 3 + 54 + 1] = MatLinCryptPrim[i * 3 + 54 + 1];
		MatLinDecrypt[PermInv[i] * 3 + 54 + 2] = MatLinCryptPrim[i * 3 + 54 + 2];
	}

	for (i = 0; i < 54; i++)
		MatLinDecrypt[i] = MatLinCrypt[i];

	free(R);
	free(MatLinCrypt);
	free(MatLinCryptPrim);
	free(Perm);
	free(PermInv);
	WriteBmpImg(BmpImg3, MatLinDecrypt);
	fclose(fileseed);

}


void colorare(unsigned char **MatImag, int x, int y, unsigned char *culoare, int HS, int WS)
{
	int i, j;
	for (i = x; i <= x + HS; i++)
	{
		MatImag[i][y] = culoare[0];
		MatImag[i][y + 1] = culoare[1];
		MatImag[i][y + 2] = culoare[2];

		MatImag[i][y + WS * 3] = culoare[0];
		MatImag[i][y + WS * 3 + 1] = culoare[1];
		MatImag[i][y + WS * 3 + 2] = culoare[2];
	}

	for (j = y; j <= y + WS * 3; j += 3)
	{
		MatImag[x][j] = culoare[0];
		MatImag[x][j + 1] = culoare[1];
		MatImag[x][j + 2] = culoare[2];

		MatImag[x + HS][j] = culoare[0];
		MatImag[x + HS][j + 1] = culoare[1];
		MatImag[x + HS][j + 2] = culoare[2];
	}


}

void templatematching(unsigned char *Imagine, unsigned char *Sablon, double prags, unsigned char *culoare, fer *fereastra, int* contor)
{
	// transform in matrice

	///////////////////////////////////////////////////////////////////////

	int WI = Imagine[19] * 256 + Imagine[18];
	int HI = Imagine[23] * 256 + Imagine[22];

	int WS = Sablon[19] * 256 + Sablon[18];
	int HS = Sablon[23] * 256 + Sablon[22];

	unsigned char **MatImag;

	MatImag = (unsigned char**)malloc(HI * sizeof(unsigned char*));


	int i, j;

	for (i = 0; i < HI; i++)
		MatImag[i] = (unsigned char*)malloc(WI * 3 * sizeof(unsigned char));



	grayscale(Imagine);
	grayscale(Sablon);

	for (i = 0; i < HI; i++)
		for (j = 0; j < WI * 3; j++)
			MatImag[i][j] = Imagine[(i)*(WI * 3) + j + 54];

	unsigned char **MatSab;

	MatSab = (unsigned char**)malloc(HS * sizeof(unsigned char*));



	for (i = 0; i < HS; i++)
		MatSab[i] = (unsigned char*)malloc(WS * 3 * sizeof(unsigned char));

	for (i = 0; i < HS; i++)
		for (j = 0; j < WS * 3; j++)
			MatSab[i][j] = Sablon[(i)*(WS * 3) + j + 54];

	///////////////////////////////////////////////////////////////////

	// functia de corelare

	///////////////////////////////////////////////////////////////////

	int n = WS * HS;
	double Scubara = 0;

	for (i = 0; i < HS; i++)
		for (j = 0; j < WS * 3; j += 3)
			Scubara += (double)MatSab[i][j];
	Scubara /= n;


	double DevStanSab = 0;
	for (i = 0; i < HS; i++)
		for (j = 0; j < WS * 3; j += 3)
			DevStanSab += (MatSab[i][j] - Scubara)*(MatSab[i][j] - Scubara);
	DevStanSab = sqrt(DevStanSab / (n - 1));


	for (i = 0; i < HI - HS; i++)
		for (j = 0; j < WI * 3 - 3 * WS; j += 3)
		{



			double corr = 0;
			double fdeimarecubara = 0;

			int k, l;

			for (k = i; k < HS + i; k++)
				for (l = j; l < WS * 3 + j; l += 3)
					fdeimarecubara += (double)MatImag[k][l];
			fdeimarecubara = (double)fdeimarecubara / (n);


			double DevStanFdei = 0;

			for (k = i; k < HS + i; k++)
				for (l = j; l < WS * 3 + j; l += 3)
				{
					unsigned char fdeimare = MatImag[k][l];
					DevStanFdei += (fdeimare - fdeimarecubara)*(fdeimare - fdeimarecubara);
				}
			DevStanFdei = (double)sqrt(DevStanFdei / (n - 1));




			for (k = i; k < HS + i; k++)
				for (l = j; l < WS * 3 + j; l += 3)
				{
					unsigned char fdeimare = MatImag[k][l];
					unsigned char sdeij = MatSab[k - i][l - j];
					corr += (double)(fdeimare - fdeimarecubara)*(sdeij - Scubara) / (DevStanFdei*DevStanSab);

				}

			corr = (double)corr / n;


			if (corr > prags)
			{

				(*contor)++;
				fereastra[(*contor) - 1].i = i;
				fereastra[(*contor) - 1].j = j;
				fereastra[(*contor) - 1].culoare[0] = culoare[0];
				fereastra[(*contor) - 1].culoare[1] = culoare[1];
				fereastra[(*contor) - 1].culoare[2] = culoare[2];
				fereastra[(*contor) - 1].correlation = corr;
				fereastra[(*contor) - 1].ok = 1;

			}



		}
	for (i = 0; i < HI; i++)
		free(MatImag[i]);
	free(MatImag);

	for (i = 0; i < HS; i++)
		free(MatSab[i]);
	free(MatSab);

}


int cmp(const void *a, const void *b)
{
	double aprim = ((fer*)a)->correlation;
	double bprim = ((fer*)b)->correlation;

	if (aprim < bprim)
		return 1;
	else
		return -1;
	return 0;
}


void ElimNonMaxime(fer* fereastra, int contor, int WS, int HS, unsigned char *Imagine)
{  
	FILE* filepragsuprapunere = fopen("prags_suprapunere.txt", "r");
	double pragsuprapunere;
	fseek(filepragsuprapunere, sizeof(float), SEEK_SET);
	fscanf(filepragsuprapunere, "%lf", &pragsuprapunere);
	fclose(filepragsuprapunere);

	int i, j;
	for (i = 0; i < contor - 1; i++)
		for (j = i + 1; j < contor; j++)
		{


			if (fereastra[i].ok == 1 && fereastra[j].ok == 1)
			{
				double suprapunere = 0;
				int ariai = HS * WS;
				int ariaj = HS * WS;
				int ariesuprapunere = 0;

				int A = max(fereastra[i].i, fereastra[j].i);
				int B = max(fereastra[i].j, fereastra[j].j);
				int C = min(fereastra[i].i + HS, fereastra[j].i + HS);
				int D = min(fereastra[i].j + WS * 3, fereastra[j].j + WS * 3);



				ariesuprapunere = (C - A + 1)*(D - B + 1) / 3;
				if (C < A || D < B)
					ariesuprapunere = 0;
			
				if ((ariai + ariaj - ariesuprapunere) != 0)
					suprapunere = (double)(ariesuprapunere) / (ariai + ariaj - ariesuprapunere);
				if (suprapunere > pragsuprapunere)
					fereastra[j].ok = 0;
			}

		}
}


int main()
{   

	//criptare
	FILE* paths = fopen("paths.txt", "r");
	char Imaginit[100], Imagcrypt[100], Imagdecrypt[100], secretkey[100];
	fscanf(paths, "%s", Imaginit);
	fscanf(paths, "%s", Imagcrypt);
	fscanf(paths, "%s", secretkey);
	fscanf(paths, "%s", Imagdecrypt);
	Crypting(Imaginit, Imagcrypt, secretkey);
	Decrypt(Imagcrypt, Imagdecrypt, secretkey);

	// template matching
	FILE* prag = fopen("prags_suprapunere.txt", "r");
	
	double prags;
	char pathzero[100];
	char pathone[100];
	char pathtwo[100];
	char paththree[100];
	char pathfour[100];
	char pathfive[100];
	char pathsix[100];
	char pathseven[100];
	char patheight[100];
	char pathnine[100];
	char pathdetection[100];

	fscanf(paths, "%s", pathzero);
	fscanf(paths, "%s", pathone);
	fscanf(paths, "%s", pathtwo);
	fscanf(paths, "%s", paththree);
	fscanf(paths, "%s", pathfour);
	fscanf(paths, "%s", pathfive);
	fscanf(paths, "%s", pathsix);
	fscanf(paths, "%s", pathseven);
	fscanf(paths, "%s", patheight);
	fscanf(paths, "%s", pathnine);
	fscanf(paths, "%s", Imaginit);
	fscanf(paths, "%s", pathdetection);
	fscanf(prag, "%lf", &prags);

	fclose(prag);


	unsigned char *cifra0 = ReadBmpImg(pathzero);
	unsigned char *cifra1 = ReadBmpImg(pathone);
	unsigned char *cifra2 = ReadBmpImg(pathtwo);
	unsigned char *cifra3 = ReadBmpImg(paththree);
	unsigned char *cifra4 = ReadBmpImg(pathfour);
	unsigned char *cifra5 = ReadBmpImg(pathfive);
	unsigned char *cifra6 = ReadBmpImg(pathsix);
	unsigned char *cifra7 = ReadBmpImg(pathseven);
	unsigned char *cifra8 = ReadBmpImg(patheight);
	unsigned char *cifra9 = ReadBmpImg(pathnine);
	unsigned char *Imagineinit = ReadBmpImg(Imaginit);

	unsigned char culoarezero[3] = { 0,0,255 };
	unsigned char culoareone[3] = { 0,255,255 };
	unsigned char culoaretwo[3] = { 0,255,0 };
	unsigned char culoarethree[3] = { 255,255,0 };
	unsigned char culoarefour[3] = { 255,0, 255 };
	unsigned char culoarefive[3] = { 255,0,0 };
	unsigned char culoaresix[3] = { 192,192,192 };
	unsigned char culoareseven[3] = { 0,140,255 };
	unsigned char culoareeight[3] = { 128,0,128 };
	unsigned char culoarenine[3] = { 0,0,128 };

	int WI = Imagineinit[19] * 256 + Imagineinit[18];
	int HI = Imagineinit[23] * 256 + Imagineinit[22];
	int WS = cifra0[19] * 256 + cifra0[18];
	int HS = cifra0[23] * 256 + cifra0[22];
	int i, j;
	fer *fereastra;                                           // aloc memorie maxima vectorului D
	fereastra = (fer*)malloc(WI*HI*3 * sizeof(fer));
	for (i = 0; i < WI*HI; i++)
		fereastra[i].ok = 0;
	int  contor = 0;


	unsigned char **MatImagInitiala;
	MatImagInitiala = (unsigned char**)malloc(HI * sizeof(unsigned char*));
	for (i = 0; i < HI; i++)
		MatImagInitiala[i] = (unsigned char*)malloc(WI * 3 * sizeof(unsigned char));

	for (i = 0; i < HI; i++)
		for (j = 0; j < WI * 3; j++)
			MatImagInitiala[i][j] = Imagineinit[(i)*(WI * 3) + j + 54];

	templatematching(Imagineinit, cifra0, prags, culoarezero, fereastra, &contor);
	templatematching(Imagineinit, cifra1, prags, culoareone, fereastra, &contor);
	templatematching(Imagineinit, cifra2, prags, culoaretwo, fereastra, &contor);
	templatematching(Imagineinit, cifra3, prags, culoarethree, fereastra, &contor);
	templatematching(Imagineinit, cifra4, prags, culoarefour, fereastra, &contor);
	templatematching(Imagineinit, cifra5, prags, culoarefive, fereastra, &contor);
	templatematching(Imagineinit, cifra6, prags, culoaresix, fereastra, &contor);
	templatematching(Imagineinit, cifra7, prags, culoareseven, fereastra, &contor);
	templatematching(Imagineinit, cifra8, prags, culoareeight, fereastra, &contor);
	templatematching(Imagineinit, cifra9, prags, culoarenine, fereastra, &contor);

	fer *aux;                                                           // copiez vectorul D intr-un aux pentru a aloca doar contor*sizeof(fereastra[0])
	aux = (fer*)malloc(contor * sizeof(fer));

	for (i = 0; i < contor; i++)
	{
		aux[i].correlation = fereastra[i].correlation;
		aux[i].culoare[0] = fereastra[i].culoare[0];
		aux[i].culoare[1] = fereastra[i].culoare[1];
		aux[i].culoare[2] = fereastra[i].culoare[2];
		aux[i].i = fereastra[i].i;
		aux[i].j = fereastra[i].j;
		aux[i].ok = fereastra[i].ok;
	}

	free(fereastra);

	qsort(aux, contor, sizeof(aux[0]), &cmp);
	ElimNonMaxime(aux, contor, WS, HS, Imagineinit);


	for (i = 0; i < contor; i++)
		if (aux[i].ok == 1)
			colorare(MatImagInitiala, aux[i].i, aux[i].j, aux[i].culoare, HS, WS);


	for (i = 0; i < HI; i++)
		for (j = 0; j < WI * 3; j++)
			Imagineinit[(i)*(WI * 3) + j + 54] = MatImagInitiala[i][j];

	WriteBmpImg(pathdetection, Imagineinit);

	fclose(paths);

	return 0;
}