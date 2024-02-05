#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <omp.h>

using namespace std;

double Comp_Isotermica(string log, double T);

string lammps_comando();

string playmol(double &sigma, double &epsilon, double &l0);

void lammps(string cmd_lammps, string arq_lammps, string arq_playmol, double sigma, double T, double P);

int MAX_ObjFunc(int n, double *comp_isot);

double **seed(int n, double xi, double xf, double yi, double yf);

double **quad_busca(double x, double y, double lx0, double ly0, double *lb, double *ub);

double comp_t(string cmd_lammps, string file_sufix, double epsilon, double sigma, double L0);

double *otimizacao(int k, int n, double epsilon, double sigma, double L0, double *lb, double *ub, string cmd_lammps, string file_sufix);

int restart_bool();

string restart_search(string cmd_lammps, double &sigma, double &epsilon, double &l0);

int num_cores();

void salvar_estado(double sigma, double epsilon, double L0, string file_sufix, double *xy_max, double f_max, double *lb, double *ub, double x_min, double x_max, double y_min, double y_max, int k, int n, int j, double **coord);

int main()
{

	srand(time(NULL));

	double epsilon, sigma, l0;
	string lammps_cmd = lammps_comando();
	string file_sufix;


	cout<<"Critical point search for diatomic molecule"<<endl;

	if(restart_bool()){
		file_sufix = restart_search(lammps_cmd,sigma,epsilon,l0);
	}else{
		file_sufix = playmol(sigma, epsilon, l0);
		double xt=comp_t(lammps_cmd,file_sufix,epsilon,sigma,l0);
		cout<<endl<<"The Isothermal Compressibility for critical point is "<<xt<<endl;
	}



    cout<<endl<<"######### EXECUTANDO MBAR EM TORNO DO PONTO CRITICO #########"<<endl;

    int N=25;
    double T[N],P,X_t;
    string name = ".input";

    FILE *file_pc = fopen("Ponto Critico.txt","r");
    char *linha = (char*)malloc(10*sizeof(char));
    fscanf(file_pc,"%[^\n]s",linha); //CABECALHO
//    fscanf(file_pc,"%s %s %s",linha,linha,linha); //epsilon sigma L0

/* TESTE */

    fscanf(file_pc,"%s",linha);
//    epsilon = atof(linha);
    
    fscanf(file_pc,"%s",linha);
//    sigma = atof(linha);
    
    fscanf(file_pc,"%s",linha);
//    l0 = atof(linha);
    
//    file_sufix = playmol(sigma, epsilon, l0);

///////////

    fscanf(file_pc,"%s",linha);
    T[N/2] = atof(linha);
    cout<<"T_c = "<<T[N/2]<<"\t";

    fscanf(file_pc,"%s",linha);
    P = atof(linha);
    cout<<"P_c = "<<P<<"\t";

    fscanf(file_pc,"%s",linha);
    X_t = atof(linha);
    cout<<"X_t = "<<X_t<<endl;

    //fclose(file_pc);

    for(int i=0;i<N/2;i++){
        T[N/2-(i+1)] = T[N/2]-2.5*(i+1); // [K]
        T[N/2+(i+1)] = T[N/2]+2.5*(i+1); // [K]
    }

    ofstream file_inp;
    string arq = "P"+to_string(P)+"_temp.inp";

#pragma omp barrier
{
    file_inp.open(arq);
    //Escrevendo arquivo temp.inp para usar no pyMBAR
    for(int i=0;i<N;i++){
        file_inp<<to_string(T[i])+"\n";
    }
    file_inp.close();

    int cores = num_cores();

#pragma omp parallel for num_threads(cores)
    for(int i=0;i<N;i++){
        //Rodando o LAMMPS
        string lammps_arq = "P"+to_string(P)+"_T"+to_string(T[i])+"_"+file_sufix+".input";
        string playmol_arq = file_sufix+".lmp_data";
        lammps(lammps_cmd,lammps_arq,playmol_arq,sigma, T[i], P);
    }
}

    string file = file_sufix+".input";
    string comando = "python3 pyMBAR.py "+to_string(P)+" "+arq+" "+file;
    system(comando.c_str());


    size_t pos = file_sufix.find("L0-");
    string mol = (file_sufix.substr(pos)).substr((file_sufix.substr(pos)).find("_")+1);
    string cmd_mv = "mv -t "+mol+" *.lmp_data *.playmol *.log *.lammps *.input *.pdf *.txt *.inp 2>/dev/null";

    system(cmd_mv.c_str());

    cout<<"TODOS OS RESULTADOS ESTAO NA PASTA \\"+mol+"\\"<<endl;

	return 0;
}

string restart_search(string cmd_lammps, double &sigma, double &epsilon, double &l0)
{
	FILE *log = fopen("search.log","r");

	fscanf(log,"%lf %lf %lf ",&sigma,&epsilon,&l0);

	char aux[100];
	fscanf(log,"%[^\n]s ",aux);
	string file_sufix = aux;

	double *xy_max= new double[2];
    int i_max;
    double f_max, *lb=new double[2], *ub=new double[2], x_min, x_max, y_min, y_max;
	fscanf(log,"%lf %lf %lf ",&xy_max[0],&xy_max[1],&f_max);
    fscanf(log,"%lf %lf ",&lb[0],&lb[1]);
	fscanf(log,"%lf %lf ",&ub[0],&ub[1]);

	fscanf(log,"%lf %lf ",&x_min,&x_max);
	fscanf(log,"%lf %lf ",&y_min,&y_max);

	int k, n, j;
	fscanf(log,"%d %d %d ",&k,&n,&j);

	double **coord = (double**)malloc(2 * sizeof(double*));
	coord[0]=(double*)malloc(n * sizeof(double)); //T
	coord[1]=(double*)malloc(n * sizeof(double)); //P
	for(int i=0; i<n; i++){
		fscanf(log,"%lf %lf ",&coord[0][i],&coord[1][i]);
	}

	double *comp_isot = new double [n];

    int cores = num_cores();

    for(j=j;j<k;j++){

#pragma omp parallel for num_threads(cores)
        for(int i=0;i<n;i++){
            string file = file_sufix;
            string lammps_arq = "P"+to_string(coord[1][i])+"_T"+to_string(coord[0][i])+"_"+file+".input";
            string playmol_arq = file+".lmp_data";
            lammps(cmd_lammps,lammps_arq,playmol_arq,sigma,coord[0][i],coord[1][i]);

            string log_lammps = "log_lammps_"+lammps_arq;

            comp_isot[i]=Comp_Isotermica(log_lammps,coord[0][i]);

#pragma omp critical
{
            FILE *file_log = fopen("OTIMIZACAO_R2W.log","a");
            fprintf(file_log,"Tc=%15.10lf \t Pc=%15.10lf \t X_t= %lf \n",coord[0][i],coord[1][i],comp_isot[i]);
            fclose(file_log);
}
        }

        i_max = MAX_ObjFunc(n, comp_isot);

        if(comp_isot[i_max]>f_max){
            f_max = comp_isot[i_max];
            xy_max[0] = coord[0][i_max];
            xy_max[1] = coord[1][i_max];
        }

        coord = quad_busca(xy_max[0],xy_max[1],x_max-x_min,y_max-y_min,lb,ub);
        x_min = coord[0][0];
        x_max = coord[0][1];
        y_min = coord[1][0];
        y_max = coord[1][1];

        coord = seed(n,x_min,x_max,y_min,y_max);
        salvar_estado(sigma, epsilon, l0, file_sufix, xy_max, f_max, lb, ub, x_min, x_max, y_min, y_max, k, n, j, coord);
    }

    ofstream file;
    file.open("Best_Pc_Tc.txt",ios_base::app);
    file<<"epsilon="+to_string(epsilon)+"\tsigma="+to_string(sigma)+"\tL0="+to_string(l0)+"\n";
    file<<"-->>-->>-->>\tTc="+to_string(xy_max[0])+"\tPc="+to_string(xy_max[1])+"\tX_t="+to_string(f_max)+"\n\n";
    file.close();

    file.open("Ponto Critico.txt");
    file<<"# epsilon\tsigma\tL0\tT_c [K]\tP_c [atm]\n";
    file<<to_string(epsilon)+"\t"+to_string(sigma)+"\t"+to_string(l0)+"\t"+to_string(xy_max[0])+"\t"+to_string(xy_max[1])+"\t"+to_string(f_max)+"\n";
    file.close();

    cout<<endl<<"The Isothermal Compressibility for critical point is "<<f_max<<endl;

    return file_sufix;
}

int restart_bool()
{
	system("ls -l > list.log");
	ifstream list;
	list.open("list.log");
	string line;
	int restart = 0;
	while(!list.eof()){
		list >> line;
		if(line.find("search.log") != string::npos){
			int opc;
			cout<<"A recent search for the critical point was interrupted."<<endl;
			cout<<"Do you want to restart it? (0 - NO / 1 - YES) ";
			cin>>opc;

			switch(opc){
				case 1:
					cout<<"YES"<<endl<<endl<<"Restarting the search..."<<endl;
					restart = 1;
					break;
				case 0:
					cout<<"NO"<<endl<<endl<<"A new search will be initiated..."<<endl;
					cout<<"Deleting files from old search..."<<endl;
					system("rm *.lmp_data *.playmol *.log *.input 2>/dev/null");
					break;
				default:
					cout<<"Invalid option! All files will be placed on a fold called: old_search"<<endl;
					system("mkdir old_search 2>/dev/null; mv -t old_search *.lmp_data *.playmol *.log 2>/dev/null");
					cout<<endl<<endl<<"A new search will be initiated..."<<endl;
					break;
			}

			break;
		}
	}
	list.close();

	return restart;
}

string playmol(double &sigma, double &epsilon, double &l0)
{
	double e_kb, R = 0.0019872036007626;
	int num;
	string molecula, massa, atomo;

	cout<<"Molecule name: ";
	cin>>molecula;

	string mkdir = "mkdir "+molecula+" 2>/dev/null";
	system(mkdir.c_str());

	cout<<"Atom: ";
	cin>>atomo;

	cout<<"Atom mass: ";
	cin>>massa;

	cout<<"Information about the Lennard-Jones potential"<<endl;

	cout<<"Epsilon/Kb [K] = ";
	cin>>e_kb;

	cout<<"sigma [A] = ";
	cin>>sigma;

	cout<<"Bond length"<<endl;

	cout<<"L0 [A] = ";
	cin>>l0;

	cout<<endl<<"How many molecules of "<<molecula<<" on simulation box: ";
	cin>>num;

	epsilon = R*e_kb;

	string file_sufix = "epsilon-"+to_string(epsilon)+"_sigma-"+to_string(sigma)+"_L0-"+to_string(l0)+"_"+molecula;
	string ext[2] = {".playmol",".lmp_data"};
	string file_name[2] = {file_sufix+ext[0],file_sufix+ext[1]};

	ofstream file;

	ostringstream epsilon_aux;
	epsilon_aux << fixed << setprecision(10) << epsilon;
	string epsilon_playmol = epsilon_aux.str();

	file.open(file_name[0]);

    file<<"define		density as 0.8765 # g/cm^3\n";
    file<<"define 		N as 300\n";
    file<<"define 		kb as 1.987204259e-3\n";
    file<<"define charge_C  as	-0.095\n";
    file<<"define charge_H  as	0.095\n";
    file<<"atom_type	C 0.06100717075	3.60    # eps 30.70 K, sig 3.60 A\n";
    file<<"atom_type	H 0.05057434839	2.36	# eps 25.45 K, sig 2.36 A\n";
    file<<"mass		C 12.011\n";
    file<<"mass       	H 11.00784\n";		
    file<<"bond_type  	C C    	1.392		# [Angstrom]\n";
    file<<"bond_type  	C H    	1.08		# [Angstrom]\n";
    file<<"angle_type 	C C C 	120		# Degrees; Energy K = 0\n";
    file<<"angle_type 	H C C 	120		# Degrees; Energy K = 0\n";
    file<<"atom	    	C1 C $charge_C\n";
    file<<"atom	    	C2 C $charge_C\n";
    file<<"atom		C3 C $charge_C\n";
    file<<"atom		C4 C $charge_C\n";
    file<<"atom		C5 C $charge_C\n";
    file<<"atom		C6 C $charge_C\n";
    file<<"atom	    	H1 H $charge_H\n";
    file<<"atom	    	H2 H $charge_H\n";
    file<<"atom		H3 H $charge_H\n";
    file<<"atom		H4 H $charge_H\n";
    file<<"atom		H5 H $charge_H\n";
    file<<"atom		H6 H $charge_H\n";
    file<<"bond     C1 C2 C6 H1\n";
    file<<"bond		C2 C3 H2\n";
    file<<"bond		C3 C4 H3\n";
    file<<"bond		C4 C5 H4\n";
    file<<"bond		C5 C6 H5\n";
    file<<"bond		C6 H6\n";
    file<<"build\n";
    file<<"12\n";
    file<<"#\n";
    file<<"  C1        0.00000        1.40272        0.00000\n";
    file<<"  H1        0.00000        2.49029        0.00000\n";
    file<<"  C2       -1.21479        0.70136        0.00000\n";
    file<<"  H2       -2.15666        1.24515        0.00000\n";
    file<<"  C3       -1.21479       -0.70136        0.00000\n";
    file<<"  H3       -2.15666       -1.24515        0.00000\n";
    file<<"  C4        0.00000       -1.40272        0.00000\n";
    file<<"  H4        0.00000       -2.49029        0.00000\n";
    file<<"  C5        1.21479       -0.70136        0.00000\n";
    file<<"  H5        2.15666       -1.24515        0.00000\n";
    file<<"  C6        1.21479        0.70136        0.00000\n";
    file<<"  H6        2.15666        1.24515        0.00000 \n";
    file<<"box         	density {0.602214*$density} # Da/A^3\n";
    file<<"packmol     	seed 1234 retry 0.95 copy mol(C1) $N action execute\n";
    file<<"write\t lammps\t"+file_name[1]+"\n";

	file.close();

	string comando = "playmol "+file_name[0];
	system(comando.c_str());

	return file_sufix;
}

string lammps_comando()
{
	system("lmp -in > output.log 2>&1");
	ifstream log;
	log.open("output.log");
	string line;
	string comando[2] = {"lmp_stable -in","lmp -in"};
	int c = 0;
	while(!log.eof()){
		log >> line;
		if(line.find("ERRO") != string::npos){
			c = 1;
		}
	}
	log.close();

	return comando[c];
}

//####################### LAMMPS EXEC ####################################
void lammps(string cmd_lammps, string arq_lammps, string arq_playmol, double sigma, double T, double P)
{
    ofstream file;
    file.open(arq_lammps);

        //Escrevendo arquivo de input do LAMMPS
    file<<"log log_lammps_"+arq_lammps+"\n\n";

    file<<"# Two sites: aa-RIGID (Nunes, 2018)\n\n";

    file<<"# --- UNIDADES, ESTADO E PARAMETROS -------------------------------------\n\n";

    file<<"units real\n";
    file<<"atom_style full\n";
    file<<"variable T equal "+to_string(T)+" # [K]\n";
    file<<"variable P equal "+to_string(P)+" # [atm]\n";
    file<<"variable sigma equal "+to_string(sigma)+" # [A]\n";
    file<<"variable rc equal 5.0*${sigma} # [A]\n\n";

    file<<"# --- ATOMOS E CAIXA DE SIMULACAO --------------------------------------- \n\n";

    file<<"pair_style lj/cut/coul/long ${rc} \n";
    file<<"pair_modify mix arithmetic tail yes \n";
    file<<"bond_style zero\n\n";
    file<<"angle_style zero\n\n";
    file<<"dihedral_style zero nocoeff\n\n";
    file<<"kspace_style pppm 1.0e-4\n\n";
    file<<"kspace_modify gewald 1.0e-4\n\n";

    file<<"# --- CONFIGURACAO INICIAL ---------------------------------------------- \n\n";

    file<<"read_data "+arq_playmol+"\n";
    file<<"velocity all create $T 261086 rot yes dist gaussian\n\n";

    file<<"# --- EQUILIBRACAOO E DEFINICAOO DO ENSEMBLE ------------------------------ \n\n";

    file<<"timestep 1.0 # [fs]\n";

    file<<"fix NVE all rigid/nve molecule\n";
    file<<"run 5000\n";
    file<<"unfix NVE\n\n";

    file<<"fix NPT all rigid/npt  molecule temp $T $T 100.0 iso $P $P 1000.0\n";
    file<<"run 500000\n\n";

    file<<"# --- EXECUCAOO DA SIMULACAO ---------------------------------------------\n\n";

    file<<"reset_timestep 0 \n";
    file<<"variable dens equal density*1000 # [kg/m^3]\n";
    file<<"thermo_style custom step temp press v_dens vol pe ke etotal &\n";
    file<<"enthalpy evdwl etail elong ebond emol\n";
    file<<"thermo 100\n";
    file<<"run 500000\n";

    file.close();

    string comando = "xterm -e "+cmd_lammps+" "+arq_lammps;
    system(comando.c_str());

}

double Comp_Isotermica(string log, double T)
{
	ifstream file;
	file.open(log.c_str());
	int ler=0,n=0,i=0;
	string line,head="Step Temp Press v_dens";
	string tail="Loop time of";

	while(getline(file,line)){
		if(ler==0 && strstr(line.c_str(),head.c_str())){
			ler=1;
		}else if(ler==1 && strstr(line.c_str(),tail.c_str())){
			ler=0;
		}else if(ler==1){
			n+=1;
		}
	}
	file.close();

	file.open(log.c_str());

	double *v = new double [n], *v2 = new double [n], *step = new double [n];
	double sum[3]={0}, v_med, v2_med;
	while(getline(file,line)){
		if(ler==0 && strstr(line.c_str(),head.c_str())){
			ler=1;
		}else if(ler==1 && strstr(line.c_str(),tail.c_str())){
			ler=0;
		}else if(ler==1){
			stringstream ss(line);
			string teste;
			int col=0;
			while(ss>>teste){
				col+=1;
				if(col==1){
					step[i] = stod(teste);
				}else if(col==5){
					v[i] = stod(teste);
					v2[i] = v[i]*v[i];
				}
			}
			sum[0]+=step[i];
			sum[1]+=v[i]*step[i];
			sum[2]+=v2[i]*step[i];
			i+=1;
		}
	}

	v_med = sum[1]/sum[0];
	v2_med= sum[2]/sum[0];
    double kb = 1.380658e-23; // [J/K]
    double beta = 1/(kb*T); // [1/J]
	return (v2_med-v_med*v_med)/v_med *  beta; // (<v^2> - <v>^2)/(<v> kb T)
}

int MAX_ObjFunc(int n, double *comp_isot)
{
    int i_max = 0;
    double f_max = -1e3;

    for(int i=0;i<n;i++){
        if(comp_isot[i]>f_max){
            i_max = i;
            f_max = comp_isot[i];
        }
    }


    return i_max;
}

double **seed(int n, double xi, double xf, double yi, double yf)
{
    double **coord = (double**)malloc(2 * sizeof(double*));
    coord[0]=(double*)malloc(n * sizeof(double)); //X
    coord[1]=(double*)malloc(n * sizeof(double)); //Y

    for(int i=0;i<n;i++){
        coord[0][i] = xi + (double)(rand()%100)/100.0 *(xf-xi); //T
        coord[1][i] = yi + (double)(rand()%100)/100.0 *(yf-yi); //P
    }

    return coord;
}

double **quad_busca(double x, double y, double lx0, double ly0, double *lb, double *ub)
{
    double **coord = (double**)malloc(2 * sizeof(double*));
    coord[0] = (double*)malloc(2 * sizeof(double)); // xi <-> xf
    coord[1] = (double*)malloc(2 * sizeof(double)); // yi <-> yf

    //x
    coord[0][0] = x - 0.7*lx0/2.0; //x_min
    coord[0][1] = x + 0.7*lx0/2.0; //x_max

    if(coord[0][0]<lb[0]){coord[0][0]=lb[0]; coord[0][1] = lb[0] + 0.7*lx0;}
    if(coord[0][1]>ub[0]){coord[0][1]=ub[0]; coord[0][0] = ub[0] - 0.7*lx0;}

    //y
    coord[1][0] = y - 0.7*ly0/2.0; //y_min
    coord[1][1] = y + 0.7*ly0/2.0; //y_max

    if(coord[1][0]<lb[1]){coord[1][0]=lb[1]; coord[1][1] = lb[1] + 0.7*ly0;}
    if(coord[1][1]>ub[1]){coord[1][1]=ub[1]; coord[1][0] = ub[1] - 0.7*ly0;}

    return coord;
}

double comp_t(string cmd_lammps, string file_sufix, double epsilon, double sigma, double L0)
{

    double *lb= new double[2];
    double *ub= new double[2];

    cout<<"Inform the temperature and pressure range for start the search..."<<endl;
    cout<<"Temperature Range [K]: \n";
    cout<<"T_MIN [K] = ";
    cin>>lb[0]; // T minimo
    cout<<"T_MAX [K] = ";
    cin>>ub[0]; // T maximo
    cout<<"\nPressure Range [atm]: \n";
    cout<<"P_MIN [atm] = ";
    cin>>lb[1]; // P minimo
    cout<<"P_MAX [atm] = ";
    cin>>ub[1]; // P maximo
    
    //cin>>lb[1]; // P central

    //ub[0]=lb[0]+100; // T_max [K]
    //ub[1]=lb[1]+100; // P_max [atm]

    //lb[0]=lb[0]-100; // T_min [K]
    //lb[1]=lb[1]-100; // P_min [atm]

    int num_int = 20;
    int gen = 16;

    double *resultado = otimizacao(num_int,gen,epsilon,sigma,L0,lb,ub,cmd_lammps,file_sufix);

    ofstream file;
    file.open("Best_Pc_Tc.txt",ios_base::app);
    file<<"epsilon="+to_string(epsilon)+"\tsigma="+to_string(sigma)+"\tL0="+to_string(L0)+"\n";
    file<<"-->>-->>-->>\tTc="+to_string(resultado[0])+"\tPc="+to_string(resultado[1])+"\tX_t="+to_string(resultado[2])+"\n\n";
    file.close();

    file.open("Ponto Critico.txt");
    file<<"# epsilon\tsigma\tL0\tT_c [K]\tP_c [atm]\n";
    file<<to_string(epsilon)+"\t"+to_string(sigma)+"\t"+to_string(L0)+"\t"+to_string(resultado[0])+"\t"+to_string(resultado[1])+"\t"+to_string(resultado[2])+"\n";
    file.close();

    return resultado[2];
}

double *otimizacao(int k, int n, double epsilon, double sigma, double L0, double *lb, double *ub, string cmd_lammps, string file_sufix)
{
    //X - T
    //Y - P
    double x_min = lb[0], x_max = ub[0], y_min = lb[1], y_max = ub[1];

    double *xy_max= new double[2];
    int i_max;
    double f_max = -1e3;

    for(int j=0;j<k;j++){
        double **coord = seed(n,x_min,x_max,y_min,y_max);

        salvar_estado(sigma, epsilon, L0, file_sufix, xy_max, f_max, lb, ub, x_min, x_max, y_min, y_max, k, n, j, coord);

        double *comp_isot = new double [n];

        int cores = num_cores();

#pragma omp parallel for num_threads(cores)
        for(int i=0;i<n;i++){
            string file = file_sufix;
            string lammps_arq = "P"+to_string(coord[1][i])+"_T"+to_string(coord[0][i])+"_"+file+".input";
            string playmol_arq = file+".lmp_data";
            lammps(cmd_lammps,lammps_arq,playmol_arq,sigma,coord[0][i],coord[1][i]);

            string log_lammps = "log_lammps_"+lammps_arq;

            comp_isot[i]=Comp_Isotermica(log_lammps,coord[0][i]);

#pragma omp critical
{
            FILE *file_log = fopen("OTIMIZACAO_R2W.log","a");
            fprintf(file_log,"Tc=%15.10lf \t Pc=%15.10lf \t X_t= %lf \n",coord[0][i],coord[1][i],comp_isot[i]);
            fclose(file_log);
}
        }

        i_max = MAX_ObjFunc(n, comp_isot);
        if(comp_isot[i_max]>f_max){
                f_max = comp_isot[i_max];
                xy_max[0] = coord[0][i_max];
                xy_max[1] = coord[1][i_max];
        }

        coord = quad_busca(xy_max[0],xy_max[1],x_max-x_min,y_max-y_min,lb,ub);
        x_min = coord[0][0];
        x_max = coord[0][1];
        y_min = coord[1][0];
        y_max = coord[1][1];
    }

    double *results = new double [3];
    results[0] = xy_max[0]; //Tc
    results[1] = xy_max[1]; //Pc
    results[2] = f_max; // X_t (Pc,Tc)

    return results;
}

void salvar_estado(double sigma, double epsilon, double L0, string file_sufix, double *xy_max, double f_max, double *lb, double *ub, double x_min, double x_max, double y_min, double y_max, int k, int n, int j, double **coord)
{
    FILE *log = fopen("search.log","w");

	fprintf(log,"%lf\t%lf\t%lf\n",sigma,epsilon,L0);

	fprintf(log,"%s\n",file_sufix.c_str());

	fprintf(log,"%lf\t%lf\t%lf\n",xy_max[0],xy_max[1],f_max);
	fprintf(log,"%lf\t%lf\n",lb[0],lb[1]);
	fprintf(log,"%lf\t%lf\n",ub[0],ub[1]);

	fprintf(log,"%lf\t%lf\n",x_min,x_max);
	fprintf(log,"%lf\t%lf\n",y_min,y_max);

	fprintf(log,"%d\t%d\t%d\n",k,n,j);

	for(int i=0; i<n; i++){
		fprintf(log,"%lf\t%lf\n",coord[0][i],coord[1][i]);
	}

	fclose(log);
}

int num_cores()
{

        int n;

        system("lscpu > CPU.lst");

        ifstream log;
        log.open("CPU.lst");
        string line;
        while(!log.eof()){
                log >> line;
                if(line.find("CPU(s):") != string::npos){
                        log >> line;
                        n = stoi(line);
                        break;
                }
        }
        log.close();

        system("rm CPU.lst");

        return n;
}
