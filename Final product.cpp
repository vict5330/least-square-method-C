#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <iomanip>
#define nmax 100

using namespace std;

//brugervejledning
void Intro ();

//Generelle funktioner
void UdskrivMatrix(double M[nmax][nmax], int n, int m);
void UdskrivVektor(double y[nmax], int n);
void UdskrivTotalMatrix(double TM[nmax][nmax+1], int n, int m);
void Transponer(double A[nmax][nmax], double AT[nmax][nmax], int n, int m);
void MatrixProdukt(double A[nmax][nmax], double AT[nmax][nmax], double MatrixProd[nmax][nmax], int n, int m);
void MatrixVektorProdukt(double A[nmax][nmax], double b[nmax], double VektorProd[nmax], int n, int m);
void TotalMatrix(double TM[nmax][nmax+1], double A[nmax][nmax], double b[nmax], int n, int m);

//Funktioner til datagrundlag a)
void HentMatrixogvektor(double A[nmax][nmax], double b[nmax], int &n, int &m, string fn);

//Funktioner til datagrundlag b)
void IndtastMatrixogvektor(double A[nmax][nmax], double b[nmax], int &n, int &m);
void GemMatrixogvektor(double A[nmax][nmax], double b[nmax], double Total[nmax][nmax+1], int n, int m, string fn);

//Funktioner til datagrundlag c)
void HentDatac(double &x1, double &deltax, double &y1, double &deltay, int &mx, int &ny, double &nc);
void Danr(double r[nmax], double nc);
void Tastabcsigma(double &ac, double &bc, double &cc, double &sigmac);
void YdreProd(double M[nmax][nmax], double a[nmax], double b[nmax], int n, int m);
void StakCols(double M[nmax][nmax], double v[nmax], int n, int m);
void XogY(double X[nmax][nmax], double Y[nmax][nmax], double x1, double deltax, double y1, double deltay, int mx, int ny);
void DanAc(double A[nmax][nmax], double x[nmax], double y[nmax], int nc);
void Danbc(double A[nmax][nmax], double b[nmax], double r[nmax], double ac, double bc, double cc, double sigmac, int nc);
//Funktioner til fordelingskontrol
void NewtonRaphsonFordelingskontrol();
double FixByTaylor(double x);
double f(double x);
void Inddelingafr(double r[nmax], int nc);

//Funktioner til datagrundlag d)
void HentDatad(double An[nmax][nmax], double A[nmax][nmax], double R[nmax][nmax], double Q[nmax][nmax], int &n, int &m, string fn);
void Danbpeogb(double An[nmax][nmax], double A[nmax][nmax], double b[nmax], int n, int m);
/*void GemMatrixogvektor(double A[nmax][nmax], double b[nmax], double Total[nmax][nmax+1], int n, int m, string fn)*/

//Funktioner til 1)
void Gauss(double TM[nmax][nmax+1], int n, int m, int &bs);
void Delvispivotering(double TM[nmax][nmax+1], int j, int n);
void BackwardsSubstitution(double TM[nmax][nmax+1], double x[nmax], int m);

//Funktioner til 2)
void GramSchimdt(double A[nmax][nmax], double Q[nmax][nmax], double R[nmax][nmax], int n, int m);
void QTbVektor(double Q[nmax][nmax], double b[nmax], double QTb[nmax], int n, int m);
/*void BackwardsSubstitution(double TM[nmax][nmax+1], double x[nmax], int m);*/

//Funktioner til 3)
double PrikProdukt(double v[nmax], double vT[nmax], int n);
double Norm(double v[nmax], int n);
void MinusVektor(double v[nmax], double vT[nmax], double Res[nmax], int n);
void HentData3(double x0[nmax], double &epsilon3, double &epsilonGolden, double &d, int &N, int m);
double f2(double A[nmax][nmax], double b[nmax], double x0[nmax], int n, int m);
void Gradient(double A[nmax][nmax], double b[nmax], double x0[nmax], double Gradientf[nmax], int n, int m);
void HesseMatrix(double A[nmax][nmax], double Hessef[nmax][nmax], int n, int m);
void InversHesseMatrix(double S[nmax][nmax], double ST[nmax][nmax], double D[nmax][nmax], double InversHessef[nmax][nmax], int n, int m);
void NewtonDirection(double InversHessef[nmax][nmax], double Gradientf[nmax], double sk[nmax], int n, int m);
double falpha(double A[nmax][nmax], double b[nmax], double alpha, double x0[nmax], double sk[nmax], int n, int m);
void Indkredsning(double A[nmax][nmax], double b[nmax], double x0[nmax], double sk[nmax], double d, double &astart, double &bstart, int n, int m);
double GoldenSection(double A[nmax][nmax], double b[nmax], double x0[nmax], double sk[nmax], double epsilonGolden, double astart, double bstart, int n, int m);
void Minimeringaff(double A[nmax][nmax], double b[nmax], double InversHessef[nmax][nmax], double Gradientf[nmax],
double v[nmax], double sk[nmax], double x0[nmax], double d, double epsilon3, double epsilonGolden, int N, int n, int m);

//Funktioner til 4)
void Tastz0epsilon4N(double z0[nmax], double &epsilon4, int &N, int m);
void PotensMetoden(double MatrixProd[nmax][nmax], double S[nmax][nmax], double D[nmax][nmax],
double ST[nmax][nmax], double SDST[nmax][nmax], double z0[nmax], double epsilon4, int n, int m, int N);

int main()
{
    char Valgabcd, ValgGemB, ValgGemD, ValgKontrolr, Valg124, Valg3, ValgIgen;
    //Generelle variabler
    double TM[nmax][nmax+1], MatrixProd[nmax][nmax], VektorProd[nmax], A[nmax][nmax], AT[nmax][nmax], b[nmax];
    int n = 0, m = 0, N, bs;

    //Variabler til datagrundlag a)
    /*double A[nmax][nmax], b[nmax]*/

    //Variabler til datagrundlag b)
    /*double A[nmax][nmax], b[nmax], TM[nmax][nmax+1]*/

    //Variabler til datagrundlag c)
    double /*A[nmax][nmax], b[nmax],*/ X[nmax][nmax], Y[nmax][nmax], r[nmax], x[nmax],
	y[nmax], ac, bc, cc, sigmac, x1, deltax, y1, deltay, nc;
    int mx, ny;

    //Variabler til datagrundlag d)
    double /*A[nmax][nmax], b[nmax],*/ An[nmax][nmax], R[nmax][nmax], Q[nmax][nmax];

    //Variabler til 1)
    /*double TM[nmax][nmax+1], x[nmax]*/

    //Variabler til 2)
    double /*A[nmax][nmax], Q[nmax][nmax], R[nmax][nmax], b[nmax], x[nmax]*/ RQTb[nmax][nmax+1], QTb[nmax];

    //Variabler til 3)
    double /*A[nmax][nmax], S[nmax][nmax], D[nmax][nmax], ST[nmax][nmax], b[nmax], */ Hessef[nmax][nmax],
	InversHessef[nmax][nmax], x0[nmax], v[nmax], Gradientf[nmax], sk[nmax], epsilon3, epsilonGolden, d;

    //Variabler til 4)
    double /*MatrixProd[nmax][nmax],*/ S[nmax][nmax], D[nmax][nmax], ST[nmax][nmax], SDST[nmax][nmax], z0[nmax], epsilon4;

    Intro ();

    do {
        cout << "\nVælg a, b, c eller d: ";
        cin >> Valgabcd;

        switch (Valgabcd)
        {
            case 'a':
            {
                cout << "\nDu har valgt at indlæse en matrix A(nxm) og en vektor b(nx1) fra en datafil." << endl;

                cout << "\nInddatafilen indlæses." << endl;

                HentMatrixogvektor(A, b, n, m, "Ab.txt");
            }
                break;

            case 'b':
            {
                cout << "\nDu har valgt at indtaste elementerne i en matrix A(nxm) og en vektor b(nx1) samt n og m." << endl;

                cout << "Indtastning af A, b, n og m." << endl;

                IndtastMatrixogvektor(A, b, n, m);

                cout << "Ønsker du at udskrive din indtastede data til en datafil?\nTast"
                " 'j' for at få udskrevet din data til en datafil og 'n' for ikke at få udskrevet din data: ";
                cin >> ValgGemB;

                if (ValgGemB == 'j')
                {
                    cout << "Du har valgt at få udskrevet din data til en datafil." << endl;

                    cout << "Udskriver n, m, A og b til en datafil." << endl;

                    GemMatrixogvektor(A, b, TM, n, m, "Abb.txt");
                }
            }
                break;

            case 'c':
            {
                cout << "\nDu har valgt at frembringe en matrix A(nxm) og en vektor b(nx1) ud fra et simuleret observationsset." << endl;

                cout << "Indtastning af x1, deltax og mx samt y1, deltay og ny. n=mx*ny dannes og r simuleres." << endl;

                HentDatac(x1, deltax, y1, deltay, mx, ny, nc);

                Danr(r, nc);

                cout << "Ønskes en fordelingskontrol af r?\nTast 'j'"
                " for at at få foretaget en fordelingskontrol af r og tast 'n' for ikke at få foretaget fordelingskontrol: ";
                cin >> ValgKontrolr;

                if (ValgKontrolr == 'j')
                {
                    cout << "Du har valgt at få foretaget en fordelingskontrol af vektoren r." << endl;

                    NewtonRaphsonFordelingskontrol();

                    Inddelingafr(r, nc);
                }

                cout << "Indtastning a, b, c og en standardafvigelse sigma hvorefter A og b dannes." << endl;

                Tastabcsigma(ac, bc, cc, sigmac);

                XogY(X, Y, x1, deltax, y1, deltay, mx, ny);

                StakCols(X, x, ny, mx);

                StakCols(Y, y, ny, mx);

                DanAc(A, x, y, nc);

                Danbc(A, b, r, ac, bc, cc, sigmac, nc);

                n = nc;

                m = 3;
            }
                break;

            case 'd':
            {
                cout << "\nDu har valgt at fremstille en matrix A(nxm) og en vektor b(nx1) ved konstruktion af overbestemte ligningssystemer Ax=b." << endl;

                cout << "Indlæsning af A8 og vælg n = 4 eller n = 8 samt indtastning af m < n og den øvre trekantsmatrix R. Q dannes og A beregnes." << endl;

                HentDatad(An, A, R, Q, n, m, "A8.txt");

                cout << "Indtastning af x* og t*. bp, e og b dannes." << endl;

                Danbpeogb(An, A, b, n, m);

                cout << "Ønsker du at udskrive dataen n, m, A og b til en datafil? "
                "Tast 'j' for at få udskrevet din data til en datafil og 'n' for ikke at få udskrevet din data: ";
                cin >> ValgGemD;

                if (ValgGemD == 'j')
                {
                    cout << "Du har valgt at få udskrevet din data til en datafil." << endl;

                    GemMatrixogvektor(A, b, TM, n, m, "Abd.txt");
                }
            }
                break;
        }

        cout << "Beregning af A^T*A og A^T*b." << endl;

        Transponer(A, AT, n, m);

        MatrixProdukt(AT, A, MatrixProd, m, n);

        cout << "\nMatrixproduktet A^T*A er givet ved:" << endl;

        UdskrivMatrix(MatrixProd, m, m);

        MatrixVektorProdukt(AT, b, VektorProd, m, n);

        cout << "Matrixvektorproduktet A^T*b er givet ved:" << endl;

        UdskrivVektor(VektorProd, m);

        cout << "Der kan vælges mellem tre muligheder til bestemmelse af x*." << endl;

        cout << "Metode 1 bestemmer x* ved løsning af normalligningerne." << endl;

        cout << "Metode 2 bestemmer x* ved QR-faktorisering." << endl;

        cout << "Metode 4 bestemmer egenløsingerne ved hjælp af potensmetoden." << endl;

        cout << "\nVælg hvilken metode du vil benytte: ";
        cin >> Valg124;

        switch (Valg124)
        {
            case '1':
            {
                cout << "\nx* findes ved løsning af normallignengerne A^T*A*x=A^T*b." << endl;

                TotalMatrix(TM, MatrixProd, VektorProd, n, m);

                Gauss(TM, n, m, bs);

                BackwardsSubstitution(TM, x, m);

                cout << "\nx* er givet ved vektoren:" << endl;

                UdskrivVektor(x, m);
            }
                break;

            case '2':
            {
                cout << "\nx* findes ved at danne A=QR og løse R*x=Q^T*b." << endl;

                GramSchimdt(A, Q, R, n, m);

                QTbVektor(Q, b, QTb, n, m);
                cout<<"Matrix Q:"<<endl;
                UdskrivMatrix(Q,n,m);
                cout<<"Matrix R:"<<endl;
                UdskrivMatrix(R,m,m);
                TotalMatrix(RQTb, R, QTb, n, m);

                BackwardsSubstitution(RQTb, x, m);

                cout << "\nx* er givet ved vektoren:" << endl;

                UdskrivVektor(x, m);
            }
                break;

            case '4':
            {
                cout << "\nEgenløsnignerne beregnes for A^T*A(mxm)." << endl;

                Tastz0epsilon4N(z0, epsilon4, N, m);

                PotensMetoden(MatrixProd, S, D, ST, SDST, z0, epsilon4, n, m, N);

                cout << "Ønsker du at vælge en metode 3 som finder x* ved at minimere f(x) ved en iterativ algoritme." << endl;

                cout << "Tast 'j' hvis du ønsker at benytte metode 3 og tast 'n' hvis ikke du ønsker at benytte metode 3: ";
                cin >> Valg3;

                if (Valg3 == 'j')
                {
                    cout << "x* findes ved minimering af f(x)." << endl;

                    HentData3(x0, epsilon3, epsilonGolden, d, N, m);

                    Gradient(A, b, x0, Gradientf, n, m);

                    HesseMatrix(A, Hessef, n, m);

                    InversHesseMatrix(S, ST, D, InversHessef, n, m);

                    NewtonDirection(InversHessef, Gradientf, sk, n, m);

                    Minimeringaff(A, b, InversHessef, Gradientf, v, sk, x0, d, epsilon3, epsilonGolden, N, n, m);
                }
            }
                break;
        }

        cout << "Ønsker du at køre programmet igen tast 'j', og hvis ikke du ønsker at køre programmet igen tast 'n': ";
        cin >> ValgIgen;

    } while (ValgIgen == 'j');

    cout << "Programmet afsluttes." << endl;

    return 0;
}

//intro
void Intro ()
{
	cout << "Velkommen til dette program, som indeholder 4 metoder til frembringelse af datagrundlag," << endl;
	cout << "samt 4 forskellige beregningsmetoder, hvilke præsenteres nedenfor. " << endl;

	cout << "\nDe fire metoder til frembringelse af datagrundlag er:" << endl;
	cout << "a) n, m og elementerne i matrix A og vektor b indlæses fra en datafil med prædefineret struktur." << endl;
	cout << "b) Brugeren indtaster n og m samt elementerne i en matrix A og en vektor b" << endl;
	cout << "c) En matrix A og en vektor b frembringes ud fra et simuleret observtionssæt (xi, yi, zi) for" << endl;
	cout << "   e lineær multipel regressionsmodel: zi = axi + byi + c + ei, med i = 1, 2, ..., n. " << endl;
	cout << "   Koefficienterne a, b og c er brugervalgte, mens ei danner en ikke-kontrollerbar tilfældig afvigelse," << endl;
	cout << "   i form af simulerede, computergenerede observationer fra en normalfordeling." << endl;
	cout << "d) EN matrix A og en vektor b fremstilles ved konstruktion af overbestemte ligningssystemer Ax=b," << endl;
	cout << "   hvor mindste kvadraters løsningen x* og mindste kvadraters fejlen: e=b-Ax er brugervalgte." << endl;

	cout << "\nHerunder ses de fire forskellige valg til beregningsmetoder:" << endl;
	cout << "1) Bestemmelse af mindste kvadraters løsningen x* findes som løsningen til normalligningerne:" << endl;
	cout << "   Ax=b  =>  A^TAx=A^Tb. Disse løses ved Gauss elimination, delvis pivotering og backwards substitution." << endl;
	cout << "2) Bestemmelse af x* sker ved løsning af: Ax=b => Rx=Q^Tb, hvor QR-faktoriseringen A=QR frembringes med" << endl;
	cout << "   modificeret Gram-Schmidt ortogonalisering af A. Ligningen Rx=Q^Tb løses ved backwards substitution og giver x*." << endl;
	cout << "3) f(x) minimeres, hvor f(x)=|b-Ax|^2 = (b-Ax)^T(b-Ax) er en kvadratisk form, og således en funktion af m variable," << endl;
	cout << "   og det gælder således, at f(x*)=min f(x). Minimeringen af f(x) foretages med en iterativ algoritme." << endl;
	cout << "4) Egenløsningerne for en matrix A bestemmes, og til dette benyttes spektralfremstillingen af A(nxn)," << endl;
	cout << "   hvor v1, v2, ..., vn er ortonormale, sådan at: A=SDS^-1 = SDS^T. Til beregning af egenløsningerne benyttes" << endl;
	cout << "   potensmetoden på de n matricer P1, P2, P3,..., Pn.";
}

//Generelle funktioner
void UdskrivMatrix(double M[nmax][nmax], int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (abs(M[i][j]) < 0.001)
            {
                M[i][j] = 0;
            }

            cout << setw(10) << M[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void UdskrivVektor(double y[nmax], int n)
{
    for(int i = 0; i < n; i++)
    {
        if (abs(y[i]) < 0.001)
        {
            y[i] = 0;
        }
        cout << setw(10) << y[i] << endl;
    }
    cout << endl;
}

void UdskrivTotalMatrix(double TM[nmax][nmax+1], int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < m+1; j++)
        {
            cout <<  setw(8) << TM[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void Transponer(double A[nmax][nmax], double AT[nmax][nmax], int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            AT[j][i] = A[i][j];
        }
    }
}

void MatrixProdukt(double AT[nmax][nmax], double A[nmax][nmax], double MatrixProd[nmax][nmax], int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            MatrixProd[i][j] = 0;

            for (int k = 0; k < m; k++)
            {
                MatrixProd[i][j] += AT[i][k]*A[k][j];
            }
        }
    }
}

void MatrixVektorProdukt(double A[nmax][nmax], double b[nmax], double VektorProd[nmax], int n, int m)
{
    for(int i = 0; i < n; i++)
    {
        VektorProd[i] = 0;

        for(int j = 0; j < m; j++)
        {
        VektorProd[i] += A[i][j]*b[j];
        }
    }
}

void TotalMatrix(double TM[nmax][nmax+1], double A[nmax][nmax], double b[nmax], int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            TM[i][j] = A[i][j];

            for (int k = 0; k <= 1; k++)
            {
                TM[i][j+1] = b[i];
            }
        }
    }
}

//Funktioner til datagrundlag a)
void HentMatrixogvektor(double A[nmax][nmax], double b[nmax], int &n, int &m, string fn)
{
    ifstream fil;

    fil.open(fn);

    fil >> n;

    fil >> m;

    for(int i = 0; i < n+1; i++)
    {
        for(int j = 0; j < m; j++)
        {
            fil >> A[i][j];
        }
        fil >> b[i];
    }
    cout << "A er givet ved: " << endl;

    UdskrivMatrix(A, n, m);

    cout << "b er givet ved: " << endl;

    UdskrivVektor(b, n);
}

//Funktioner til datagrundlag b)
void IndtastMatrixogvektor(double A[nmax][nmax], double b[nmax], int &n, int &m)
{
    cout << "\nDu skal nu indtaste en vilkårlig matrix og og en vilkårlig vektor" << endl;

    cout << "\nDu bedes indtaste n-rækker og m-søjler for nxm-matricen" << endl;

    cout << "\nn = ";
    cin >> n;

    cout << "\nm = ";
    cin >> m;

    cout << "\nMatricen 'A' indtastes nu" << endl;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            cout<< "Indtast det " << i+1 << "," << j+1 << " element i matricen: ";
            cin >> A[i][j];
        }
    }

    cout << "\nVektoren 'b' indtastes nu" << endl;

    for(int i = 0; i < n; i++)
    {
        cout << "Indtast det " << i+1 << " element i vektoren: ";
        cin >> b[i];
    }

    cout << "A er givet ved: " << endl;

    UdskrivMatrix(A, n, m);

    cout << "b er givet ved: " << endl;

    UdskrivVektor(b, n);
}

void GemMatrixogvektor(double A[nmax][nmax], double b[nmax], double TM[nmax][nmax+1], int n, int m, string fn)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            TM[i][j] = A[i][j];

            for (int k = 0; k <= 1; k++)
            {
                TM[i][j+1] = b[i];
            }
        }
    }

    ofstream fil;

    fil.open(fn);

    fil << n << endl;

    fil << m << endl;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            if(j > 0)
            {
                fil << " ";
            }
            fil << TM[i][j];
        }
        fil << endl;
    }
    fil.close();
}

//Funktioner til datagrundlag c)
void HentDatac(double &x1, double &deltax, double &y1, double &deltay, int &mx, int &ny, double &nc)
{
    cout << "Du bedes først indtaste x1, delta-x og mx." << endl;

    cout << "x1: ";
    cin >> x1;

    cout << "delta-x: ";
    cin >> deltax;

    cout << "mx: ";
    cin >> mx;

    cout << "\nDu bedes nu indtaste y1, delta-y og ny." << endl;

    cout << "y1: ";
    cin >> y1;

    cout << "delta-y: ";
    cin >> deltay;

    cout << "ny: ";
    cin >> ny;

    cout << "\nBeregner n som er ny*mx:" << endl;

    nc = ny*mx;

    cout << "n = " << nc << endl;

    cout << "\nNu simuleres r-vektoren:" << endl;
}

void Danr(double r[nmax], double nc)
{
    double ra[nmax], rb[nmax], xa[nmax], xb[nmax];
    int l = 1, k = 1;

    srand(time(0));

    for(int i = 0; i < nc/2; i++)
    {
        xa[i] = (rand()%10001)/10000.0;
        xb[i] = (rand()%10001)/10000.0;
    }

    for(int i = 0; i < nc/2; i++)
    {
        ra[i] = sqrt((-2*log(xa[i])))*cos(2*M_PI*xb[i]);
        rb[i] = sqrt((-2*log(xa[i])))*sin(2*M_PI*xb[i]);
    }

    for(int i = 0; i <= nc; i = i+2)
    {
        for(int j = 0; j < k; j++)
        {
            r[i] = ra[j];
        }
        k++;
    }

    for(int i = 1; i <= nc; i = i+2)
    {
        for(int j = 0; j < l; j++)
        {
            r[i] = rb[j];
        }
        l++;
    }

    for(int i = 0; i < nc; i++)
    {
        cout << "r[" << i+1 << "] = " << r[i] << endl;
    }
    cout << endl;
}

void Tastabcsigma(double &ac, double &bc, double &cc, double &sigmac)
{
    cout << "Der skal nu dannes en vektor ud fra koordinaterne a, b og c" << endl;

    cout << "Tast 1.koordinaten a: ";
    cin >> ac;

    cout << "Tast 2.koordinaten b: ";
    cin >> bc;

    cout << "Tast 3. koordinaten c: ";
    cin >> cc;

    cout << "Indtast også en standardafvigelse sigma større end 0: ";
    cin >> sigmac;
}

void YdreProd(double M[nmax][nmax], double a[nmax], double b[nmax], int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            M[i][j] = a[i]*b[j];
        }
    }
}

void StakCols(double M[nmax][nmax], double v[nmax], int n, int m)
{
    int a;

    for(int i = 0; i < m ; i++)
    {
        a = i*n;

        for(int j = 0; j < n ; j++)
        {
            v[a+j] = M[j][i];
        }
    }
}

void XogY(double X[nmax][nmax], double Y[nmax][nmax], double x1, double deltax, double y1, double deltay, int mx, int ny)
{
    double xT[nmax], yT[nmax], eTx[nmax], eTy[nmax];

    for(int i = 0; i < mx; i++)
    {
        xT[i] = x1+i*deltax;
    }

    for(int i = 0; i < ny; i++)
    {
        yT[i] = y1+i*deltay;
    }

    for(int i = 0; i < mx; i++)
    {
        eTx[i] = 1;
    }

    for(int i = 0; i < ny; i++)
    {
        eTy[i] = 1;
    }

    cout << "\ne^T for mx = ";
    for(int i = 0; i < mx; i++)
    {
        cout << eTx[i] << " ";
    }

    cout << "\ne^T for ny = ";
    for(int i = 0; i < ny; i++)
    {
        cout << eTy[i] << " ";
    }

    cout << "\nx^T = ";
    for(int i = 0; i < mx; i++)
    {
        cout << xT[i] << " ";
    }

    cout << "\ny^T = ";
    for(int i = 0; i < ny; i++)
    {
        cout << yT[i] << " ";
    }
    cout << endl;

    YdreProd(X,eTy,xT,ny,mx);

    cout << "X-matricen er givet ved:" << endl;

    UdskrivMatrix(X, ny, mx);

    YdreProd(Y, yT, eTx, ny, mx);

    cout << "Y-matricen er givet ved:" << endl;

    UdskrivMatrix(Y, ny, mx);
}

void DanAc(double A[nmax][nmax], double x[nmax], double y[nmax], int nc)
{
    for (int i = 0; i < nc; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A[i][0] = x[i];
            A[i][1] = y[i];
            A[i][2] = 1;
        }
    }

    cout << "\nMatricen A er givet ved:" << endl;

    UdskrivMatrix(A, nc, 3);
}

void Danbc(double A[nmax][nmax], double b[nmax], double r[nmax], double ac, double bc, double cc, double sigmac, int nc)
{
    double epsilonc[nmax], z0[nmax];

    for (int i = 0; i < nc; i++)
    {
        z0[i] = ac*A[i][0]+bc*A[i][1]+cc*A[i][2];
    }

    for (int i = 0; i < nc; i++)
    {
        epsilonc[i] = sigmac*r[i];
    }

    for (int i = 0; i < nc; i++)
    {
        b[i] = z0[i]+epsilonc[i];
    }

    cout << "\nVektoren b er givet ved:" << endl;

    UdskrivVektor(b, nc);
}

//Funktioner til fordelingskontrol
void NewtonRaphsonFordelingskontrol()
{
    double c, xny, xgl, epsilon = 10e-8, x0[7] = {-1.2, -0.7, -0.4, 0, 0.4, 0.7, 1.2};
    int j = 0, k = 0, i = 0, N =1000;

    cout << "\nNewton-Raphson udført på Phi(xk) = 0.125*k for iterationer af k fra 1 til 7:" << endl;

    do
    {
        k++;
        c = 0.125*k;

        xny = x0[i];

        do
        {
            j++;

            xgl = xny;

            xny = xgl-(FixByTaylor(xgl)-c)/f(xgl);

        } while ((fabs(FixByTaylor(xny)-c) > epsilon) && (j < N));

        if (fabs(FixByTaylor(xny)-c) > epsilon)
        {
            cout << "Der kunne ikke findes et nulpunkt." << endl;
        }
        else
        {
            cout << "Algortimen konvergerede i nulpunktet x" << k << "= " << xny << endl;

            cout << "Phi(" << xny << ") = " << FixByTaylor(xny) << endl;
        }

        i++;

    } while (k < 7);
}

double FixByTaylor(double x)
{
    double sum = x, term = x, fac = 0, f;
    int k = 1;

    do
    {
        f = sum;

        k = k + 2;

        fac = fac + 1;

        term = term * (-0.5) * ((x*x)/k) * (k-2) * (1/fac);

        sum = sum + term;
    }
    while(f != sum);

    sum = 1/sqrt(2*M_PI)*sum+0.5;

    return sum;
}

double f(double x)
{
    return (1/sqrt(2*M_PI))*exp((-1/2)*pow(x, 2));
}

void Inddelingafr(double r[nmax], int nc)
{
    int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0, i7 = 0, i8 = 0;

    for (int i = 0; i < nc; i++)
    {
        if (r[i] < -1.15035)
        {
            i1++;
        }

        if (r[i] > -1.15035 && r[i] < -0.67449)
        {
            i2++;
        }

        if (r[i] > -0.67449 && r[i] < -0.318639)
        {
            i3++;
        }

        if (r[i] > -0.318639 && r[i] < 0)
        {
            i4++;
        }

        if (r[i] > 0 && r[i] < 0.318639)
        {
            i5++;
        }

        if (r[i] > 0.318639 && r[i] < 0.67449)
        {
            i6++;
        }

        if (r[i] > 0.67449 && r[i] < 1.15035)
        {
            i7++;
        }

        if (r[i] > 1.15035)
        {
            i8++;
        }
    }

    cout << "\nVærdierne af de forskellige r er inddelt i intervallerne således: " << endl;

    cout << "Antal i 1. interval: " << i1 << endl;

    cout << "Antal i 2. interval: " << i2 << endl;

    cout << "Antal i 3. interval: " << i3 << endl;

    cout << "Antal i 4. interval: " << i4 << endl;

    cout << "Antal i 5. interval: " << i5 << endl;

    cout << "Antal i 6. interval: " << i6 << endl;

    cout << "Antal i 7. interval: " << i7 << endl;

    cout << "Antal i 8. interval: " << i8 << endl;

    cout << "Total: " << i1+i2+i3+i4+i5+i6+i7+i8 << endl;
}

//Funktioner til datagrundlag d)
void HentDatad(double An[nmax][nmax], double A[nmax][nmax], double R[nmax][nmax], double Q[nmax][nmax], int &n, int &m, string fn)
{
    double D[nmax][nmax], v[nmax], PrikProdukt;

    ifstream fil;
    fil.open(fn);

    cout << "\nVælg n lig 4 eller n lig 8: ";
    cin >> n;

    cout << "\nIndtast m, hvor m < n: ";
    cin >> m;

        for(int i = 0; i < 8; i++)
        {
            for(int j = 0; j < 8; j++)
            {
                fil >> An[i][j];
            }
        }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            v[j] = An[j][i];
        }

        PrikProdukt = 0;

        for (int j = 0; j < n; j++)
        {
            PrikProdukt += v[j]*v[j];
        }

        D[i][i] = 0;

        D[i][i] = 1/sqrt(PrikProdukt);
    }

    MatrixProdukt(An, D, Q, n, m);

    cout << "\nIndtast den øvre trekants matrix R: " << endl;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << "Indtast det " << i+1 << "," << j+1 << ". element i R: ";
            cin >> R[i][j];
        }
    }

    MatrixProdukt(Q, R, A, n, m);

    cout << "A" << n << " er givet ved:" << endl;

    UdskrivMatrix(An, n, n);

    cout << "Den ortonormale matrix Q er givet ved:" << endl;

    UdskrivMatrix(Q, n, m);

    cout << "Den øvre trekantsmatrix R er givet ved:" << endl;

    UdskrivMatrix(R, m, m);

    cout << "Matricen A er givet ved:" << endl;

    UdskrivMatrix(A, n, m);
}

void Danbpeogb(double An[nmax][nmax], double A[nmax][nmax], double b[nmax], int n, int m)
{
    double x[nmax], t[nmax], bp[nmax], e[nmax];

    cout << "Indtast x*-vektoren." << endl;

    for (int i = 0; i < m; i++)
    {
        cout << "Indtast det " << i+1 << ". element i x*: ";
        cin >> x[i];
    }

    cout << "Indtast t*-vektoren." << endl;

    for (int i = 0; i < n-m; i++)
    {
        cout << "Indtast det " << i+1 << ". element i t*: ";
        cin >> t[i];
    }

    MatrixVektorProdukt(A, x, bp, n, m);

    for (int i = 0; i < n; i++)
    {
        e[i] = 0;
        for (int j = 0; j < m; j++)
        {
            e[i] += t[j]*An[i][j+m];
        }
    }

    for (int i = 0; i < n; i++)
    {
        b[i] = bp[i]+e[i];
    }

    cout << "Vektoren x* er givet ved:" << endl;

    UdskrivVektor(x, m);

    cout << "Vektoren t* er givet ved:" << endl;

    UdskrivVektor(t, n-m);

    cout << "Vektoren bp er givet ved:" << endl;

    UdskrivVektor(bp, n);

    cout << "Vektoren e er givet ved:" << endl;

    UdskrivVektor(e, n);

    cout << "Vektoren b er givet ved:" << endl;

    UdskrivVektor(b, n);
}

//Funktioner til 1)
void Gauss (double TM[nmax][nmax+1], int n, int m, int &bs)
{
    double factor, epsilon;

    epsilon = 10e-8;

    bs = 1;

    for(int i = 0; i <= n-2; i++)
    {
        Delvispivotering (TM, i, n);

        if(fabs(TM[i][i]) < epsilon)
        {
            bs = 0;

            break;
        }

        for(int j = i+1; j <= n; j++)
        {
            factor = -TM[j][i]/TM[i][i];

            TM[j][i]=0;

            for(int k = i+1; k <= n; k++)
            {
                TM[j][k]=TM[j][k]+factor*TM[i][k];
            }
        }

        if (fabs(TM[n-1][n-1]) < epsilon)
        {
            bs=0;
        }
    }
}

void Delvispivotering(double TM[nmax][nmax+1], int j, int n)
{
    int maxr = j;

    for(int i = j+1; i < n; i++)
    {
        if(fabs(TM[maxr][j]) < fabs(TM[i][j]))
        {
        maxr = i;
        }
    }
    for(int k = j; k <= n; k++)
    {
        double temp = TM[j][k];

        TM[j][k] = TM[maxr][k];

        TM[maxr][k] = temp;
    }
}

void BackwardsSubstitution(double TM[nmax][nmax+1], double x[nmax], int m)
{
    double sum;

    x[m-1] = TM[m-1][m]/TM[m-1][m-1];

    for(int i = m-2; i >= 0; i--)
    {
        sum = 0;

        for(int j = i+1; j <= m-1; j++)
        {
            sum = sum+TM[i][j]*x[j];
        }
        x[i] = (TM[i][m]-sum)/TM[i][i];
    }
}

//Funktioner til 2)
void GramSchimdt(double A[nmax][nmax], double Q[nmax][nmax], double R[nmax][nmax], int n, int m)
{
    for (int i = 0; i < m; i++)
    {
        R[i][i] = 0;

        for (int j = 0; j < n; j++)
        {
        R[i][i] += A[j][i]*A[j][i];
        }

        R[i][i] = sqrt(R[i][i]);

        for (int j = 0; j < n; j++)
        {
            Q[j][i] = A[j][i]/R[i][i];
        }

        for(int k = i+1; k < n; k++)
        {
            R[i][k] = 0;

            for(int j = 0; j < n; j++)
            {
                R[i][k] += Q[j][i]*A[j][k];
            }

            for (int j = 0; j < n; j++)
            {
                A[j][k] -= R[i][k]*Q[j][i];
            }
        }
    }
}

void QTbVektor(double Q[nmax][nmax], double b[nmax], double QTb[nmax], int n, int m)
{
    double QT[nmax][nmax];

    Transponer(Q, QT, n, m);

    for (int i = 0; i < m; i++)
    {
        QTb[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            QTb[j] += QT[j][i]*b[i];
        }
    }
}

//Funktioner til 3)
double PrikProdukt(double v[nmax], double vT[nmax], int n)
{
    double sum = 0;

    for (int i = 0; i < n; i++)
    {
        sum += v[nmax]*vT[nmax];
    }

    return sum;
}

double Norm(double v[nmax], int n)
{
    return sqrt(PrikProdukt(v, v, n));
}

void MinusVektor(double v[nmax], double vT[nmax], double Res[nmax], int n)
{
    for (int i = 0; i < n; i++)
    {
        Res[i] = v[i]-vT[i];
    }
}

void HentData3(double x0[nmax], double &epsilon3, double &epsilonGolden, double &d, int &N, int m)
{
    cout << "Indtast startvektoren x0." << endl;

    for (int i = 0; i < m; i++)
    {
        cout << "Indtast det " << i+1 << ". element i x0: ";
        cin >> x0[i];
    }

    cout << "Indtast tolerancen for stopkriteriet i algoritmen: ";
    cin >> epsilon3;

    cout << "Indtast det maksimale antal iterationer i algoritmen: ";
    cin >> N;

    cout << "Indtast startværdien for indkredsningsalgoritmen, som finder et startinterval for Golden Section algoritmen større end 0: ";
    cin >> d;

    cout << "Indtast tolerancen i stopkriteriet for Golden Section algoritmen: ";
    cin >> epsilonGolden;
}

double f2(double A[nmax][nmax], double b[nmax], double x0[nmax], int n, int m)
{
    double Ax[nmax], Funktionsværdi = 0;

    MatrixVektorProdukt(A, x0, Ax, n, m);

    for (int i = 0; i < n; i++)
    {
        Funktionsværdi += (b[i]-Ax[i])*(b[i]-Ax[i]);
    }

    return Funktionsværdi;
}

void Gradient(double A[nmax][nmax], double b[nmax], double x0[nmax], double Gradientf[nmax], int n, int m)
{
    double Ax0[nmax], AT[nmax][nmax], ATAx0[nmax], ATb[nmax], Led1[nmax], Led2[nmax];

    MatrixVektorProdukt(A, x0, Ax0, n, m);

    Transponer(A, AT, n, m);

    MatrixVektorProdukt(AT, Ax0, ATAx0, m, n);

    for (int i = 0; i < m; i++)
    {
        Led1[i] = 2*ATAx0[i];
    }

    MatrixVektorProdukt(AT, b, ATb, m, n);

    for (int i = 0; i < m; i++)
    {
        Led2[i] = 2*ATb[i];
    }

    for (int i = 0; i < m; i++)
    {
        Gradientf[i] = Led1[i]-Led2[i];
    }
}

void HesseMatrix(double A[nmax][nmax], double Hessef[nmax][nmax], int n, int m)
{
    double AT[nmax][nmax], ATA[nmax][nmax];

    Transponer(A, AT, n, m);

    MatrixProdukt(AT, A, ATA, n, m);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Hessef[i][j] = 2*ATA[i][j];
        }
    }
}

void InversHesseMatrix(double S[nmax][nmax], double ST[nmax][nmax], double D[nmax][nmax], double InversHessef[nmax][nmax], int n, int m)
{
    double DInvers[nmax][nmax], DInversS[nmax][nmax], STDInversS[nmax][nmax];

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (j == i)
            {
                DInvers[j][j] = 1/D[j][j];
            }
            else
            {
                DInvers[j][i] = 0;
            }
        }
    }

    MatrixProdukt(S, DInvers, DInversS, m, m);

    MatrixProdukt(DInversS, ST, STDInversS, m, m);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            InversHessef[i][j] = 0.5*STDInversS[i][j];
        }
    }
}

void NewtonDirection(double InversHessef[nmax][nmax], double Gradientf[nmax], double sk[nmax], int n, int m)
{
    double NegInversHessef[nmax][nmax], sk1[nmax], Normsk1 = 0, Normsk2 = 0;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            NegInversHessef[i][j] = -1*InversHessef[i][j];
        }
    }

    MatrixVektorProdukt(NegInversHessef, Gradientf, sk1, m, m);

    for (int i = 0; i < m; i++)
    {
        Normsk1 += sk1[i]*sk1[i];
    }

    Normsk2 = sqrt(Normsk1);

    for (int i = 0; i < m; i++)
    {
        sk[i] = sk1[i]/Normsk2;
    }
}

double falpha(double A[nmax][nmax], double b[nmax], double alpha, double x0[nmax], double sk[nmax], int n, int m)
{
    double v[nmax];

    for (int i = 0; i < n; i++)
    {
        v[i] = x0[i]+alpha*sk[i];
    }

    return f2(A, b, v, n, m);
}

void Indkredsning(double A[nmax][nmax], double b[nmax], double x0[nmax],
double sk[nmax], double d, double &astart, double &bstart, int n, int m)
{
    double alpha1 = 0, alpha2, alpha3;

    alpha2 = alpha1+d;

    while (falpha(A, b, alpha2, x0, sk, n, m) >= falpha(A, b, alpha1, x0, sk, n, m))
    {
        d = 0.1*d;

        alpha2 = alpha1+d;
    }

    d = 2*d;

    alpha3 = alpha2+d;

    while (falpha(A, b, alpha2, x0, sk, n, m) >= falpha(A, b, alpha3, x0, sk, n, m))
    {
        alpha1 = alpha2;

        alpha2 = alpha3;

        d = 2*d;

        alpha3 = alpha2+d;
    }

    astart = alpha1;

    bstart = alpha3;
}

double GoldenSection(double A[nmax][nmax], double b[nmax], double x0[nmax]
,double sk[nmax], double epsilonGolden, double astart, double bstart, int n, int m)
{
    double c, alpha1 = astart, alpha2, alpha3, alpha4 = bstart;

    c = (sqrt(5)-1)/2;

    alpha2 = alpha1+(1-c)*(alpha4-alpha1);

    alpha3 = alpha4-(1-c)*(alpha4-alpha1);

    while (abs(alpha4-alpha1) >= epsilonGolden)
    {
        if (falpha(A, b, alpha2, x0, sk, n, m) >= falpha(A, b, alpha3, x0, sk, n, m))
        {
            alpha1 = alpha2;

            alpha2 = alpha3;

            alpha3 = alpha4-(1-c)*(alpha4-alpha1);
        }
        else
        {
            alpha4 = alpha3;

            alpha3 = alpha2;

            alpha2 = alpha1+(1-c)*(alpha4-alpha1);
        }
    }

    return 0.5*(alpha4+alpha1);
}

void Minimeringaff(double A[nmax][nmax], double b[nmax], double InversHessef[nmax][nmax],
double Gradientf[nmax], double v[nmax], double sk[nmax], double x0[nmax],
double d, double epsilon3, double epsilonGolden, int N, int n, int m)
{
    double xny[nmax], xgl[nmax], Res[nmax], alpha, astart, bstart;
    int k = 0;

    for (int i = 0; i < m; i++)
    {
        xny[i] = x0[i];
    }

    do
    {
        for (int i = 0; i < m; i++)
        {
            xgl[i] = xny[i];
        }

        NewtonDirection(InversHessef, Gradientf, sk, n, m);

        Indkredsning(A, b, xgl, sk, d, astart, bstart, n, m);

        alpha = GoldenSection(A, b, xgl, sk, epsilonGolden, astart, bstart, n, m);

        for (int i = 0; i < m; i++)
        {
            xny[i] = xgl[i]+alpha*sk[i];
        }

        k++;

        MinusVektor(xny, xgl, Res, m);

    } while (Norm(Res, n) >= epsilon3 && Norm(Gradientf, n) >= epsilon3 && k != N);

    cout << "x* er givet ved vektoren:" << endl;

    UdskrivVektor(xny, m);
}

//Funktioner til 4)
void Tastz0epsilon4N(double z0[nmax], double &epsilon4, int &N, int m)
{
    cout << "\nIndtast startvektoren z0 med n rækker." << endl;

    for (int i = 0; i < m; i++)
    {
        cout << "Indtast det " << i+1 << ". element i z0: ";
        cin >> z0[i];
    }

    cout << "Indtast tolerancen for stopkriteriet: ";
    cin >> epsilon4;

    cout << "Indtast den øvre grænse for antallet af tilladte iterationer: ";
    cin >> N;
}

void PotensMetoden(double MatrixProd[nmax][nmax], double S[nmax][nmax], double D[nmax][nmax],
	double ST[nmax][nmax], double SDST[nmax][nmax], double z0[nmax], double epsilon4, int n, int m, int N)
{
    double wwT[nmax][nmax], Pny[nmax][nmax], Pgl[nmax][nmax], yny[nmax], ygl[nmax],
	v[nmax], Lny, Lgl, yi = 0, normyi, y, normy, PrikProd;
    int i, k;
    bool løsning;

    for (int j = 0; j < m; j++)
    {
        for (int h = 0; h < m; h++)
        {
            Pny[j][h] = MatrixProd[j][h];
        }
    }

    for (k = 0; k < m; k++)
    {
        for (int j = 0; j < m; j++)
        {
            for (int h = 0; h < m; h++)
            {
                Pgl[j][h] = Pny[j][h];
            }
        }

        i = 0;

        Lny = z0[0];

        for (int j = 0; j < m; j++)
        {
            if (abs(z0[j]) > Lny)
            {
                Lny = z0[j];
            }
        }

        for (int j = 0; j < m; j++)
        {
            yny[j] = (1/Lny)*z0[j];
        }

        do
        {
            i++;

            Lgl = Lny;

            for (int j = 0; j < m; j++)
            {
                ygl[j] = yny[j];
            }

            MatrixVektorProdukt(Pgl, ygl, z0, m, m);

            Lny = z0[0];

            for (int j = 0; j < m; j++)
            {
                if (abs(z0[j]) > Lny)
                {
                    Lny = z0[j];
                }
            }

            for (int j = 0; j < m; j++)
            {
                yny[j] = (1/Lny)*z0[j];
            }

            for (int j = 0; j < m; j++)
            {
                yi += pow(yny[j], 2)-pow(ygl[j], 2);
            }

            normyi = sqrt(yi);

        } while (abs(Lny-Lgl) >= epsilon4 && normyi >= epsilon4 && i != N);

            if (i < N)
            {
                løsning = true;
            }

            if (løsning)
            {
                cout << "\nLøsningskriteriummet er opfyldt." << endl;

                cout << "Den numerisk største tilnærmede egenværdi efter " << i << " iterationer er givet ved: " << Lny << endl;

                y = 0;

                for (int j = 0; j < m; j++)
                {
                    y += pow(yny[j], 2);
                }

                normy = sqrt(y);

                for (int j = 0; j < m; j++)
                {
                    v[j] = yny[j]/normy;
                }

                for (int j = 0; j < m; j++)
                {
                    S[j][k] = v[j];
                }

                for (int j = 0; j < m; j++)
                {
                    if (j == k)
                    {
                        D[j][j] = Lny;
                    }
                    else
                    {
                        D[j][k] = 0;
                    }
                }

                YdreProd(wwT, yny, yny, m, m);

                PrikProd = 0;

                for (int j = 0; j < m; j++)
                {

                    PrikProd += yny[j]*yny[j];
                }

                for (int j = 0; j < m; j++)
                {
                    for (int h = 0; h < m; h++)
                    {
                        Pny[j][h] = 0;

                        Pny[j][h] = Pgl[j][h]-(Lny*(wwT[j][h]/PrikProd));
                    }
                }

                cout << "Den tilhørende tilnærmede egenvektor til egenværdien er:" << endl;

                UdskrivVektor(v, m);

                cout << "P" << k+2 << " er givet ved:" << endl;

                UdskrivMatrix(Pny, m, m);

            }
            else
            {
                cout << "\nLøsningskriteriummet er ikke opfyldt." << endl;
            }
    }
    cout << "Matricen S, der indeholder de ortonormale egenvektorer er givet ved:" << endl;

    UdskrivMatrix(S, m, m);

    cout << "Matricen D, der er en diagonalematrix og indeholder egenværdierne til matricen A^T*A, er givet ved:" << endl;

    UdskrivMatrix(D, m, m);

    MatrixProdukt(S, D, MatrixProd, m, m);

    Transponer(S, ST, m, m);

    MatrixProdukt(MatrixProd, ST, SDST, m, m);

    cout << "S*D*S^T skulle gerne give A^T*A:" << endl;

    UdskrivMatrix(SDST, m, m);
}
