#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // rand, srand
#include <time.h>       /* time */
#include <sstream>
#include <string>
#include <vector>

using namespace std;

const int LONG_CROMOSOMA = 30;
const double COEF = (1 << LONG_CROMOSOMA) - 1;

const int MAX_POP = 30;
const int MAX_GEN = 100;

const double pCross = 0.6;
const double pMutation = 0.01;


// usar escalado lineal para evitar perdida temprana de diversidad genetica?
const double ESCALAR = true;
const double fmultiple = 2.5;

// mostrar cada individuo?
const double MOSTRAR_DETALLE = false;

class Individuo{
public:
    bool cromosoma[LONG_CROMOSOMA];
    int longitud;
    double fitness;
    double scaled_fitness = -1;
    double dec;
    int idx_padre1 = 0;
    int idx_padre2 = 0;
    int pos_cross = 0;
    int cant_mutaciones;

	Individuo() {
		cant_mutaciones = 0;
	}

    double decodificar(){
        long total = 0;
        for(int i = 0; i < LONG_CROMOSOMA; i++){
            total = total * 2 + this->cromosoma[i]; // (int)True == 1 por estandar
        }
        return total / COEF;
    }


    double calcularFitness(){
        return pow(this->dec,10) ;
    }

    void Actualizar(){
        this->dec = decodificar();
        this->fitness = calcularFitness();
        this->scaled_fitness = this->fitness;
    }

};

 Individuo* GenerarRandom(){
    Individuo* ind = new Individuo();
    for(int i = 0; i < LONG_CROMOSOMA; i++){
        ind->cromosoma[i] = rand() & 1;
    }
    ind->Actualizar();
    return ind;
}

// elige un booleano aleatorio segun probabilidad p
bool flip(double p){
    return rand() < RAND_MAX * p;
}

bool mutar(bool alelo,int& cont_mut){
    bool m = flip(pMutation);
    cont_mut += m;
    return alelo ^ m;
}

float random_float(){
    return rand() / (double)RAND_MAX;
}

void Crossover(vector<Individuo*>& viejaGen,int idx_p1, int idx_p2,vector<Individuo*>& nuevaGen){
    Individuo* padre1 = viejaGen[idx_p1];
    Individuo* padre2 = viejaGen[idx_p2];
    int cross;
    if (flip(pCross)){
        // si hay crossover elije una posicion al azar donde empezar el crossover entre la 2da y la ultima
        // (con cierta probabilidad de mutación)
        cross = 1 + rand() % (LONG_CROMOSOMA - 1);
    }
    else{
        // sino no hace el crossover pasa ambos padres a la siguiente generacion con
        // solamente la probabilidad de mutación. Establesco cross en len para que
        // haga una copia entera
        cross = LONG_CROMOSOMA;
    }
    Individuo* hijo1 = new Individuo();
    Individuo* hijo2 = new Individuo();
    for(int i = 0; i < cross; i++){
        hijo1->cromosoma[i] = mutar(padre1->cromosoma[i],hijo1->cant_mutaciones);
        hijo2->cromosoma[i] = mutar(padre2->cromosoma[i],hijo2->cant_mutaciones);
    }
    for(int i = cross; i < LONG_CROMOSOMA; i++){
        hijo1->cromosoma[i] = mutar(padre2->cromosoma[i],hijo1->cant_mutaciones);
        hijo2->cromosoma[i] = mutar(padre1->cromosoma[i],hijo2->cant_mutaciones);
    }
    hijo1->Actualizar();
    hijo2->Actualizar();

    hijo1->pos_cross = cross;
    hijo1->idx_padre1 = idx_p1;
    hijo1->idx_padre2 = idx_p2;

    hijo2->pos_cross = cross;
    hijo2->idx_padre2 = idx_p1;
    hijo2->idx_padre1 = idx_p2;

    nuevaGen.push_back(hijo1);
    nuevaGen.push_back(hijo2);
};

string CromosomaString(Individuo* ind){
    std::stringstream ss;
    for(int i = 0; i < LONG_CROMOSOMA; i++){
        ss << ind->cromosoma[i];
    }
    return ss.str();
}


// retorna la posicion a elegir segun ruleta entre 0 y len - 1
int Select(const vector<Individuo*>& pop, double sum_scaled_fitness ){
    double suma = 0;
    double punto_ruleta = random_float() * sum_scaled_fitness;

    for(vector<Individuo*>::size_type i = 0; i < pop.size(); i++){
        suma += pop[i]->scaled_fitness;
        if (suma >= punto_ruleta) return i;
    }
    return pop.size() - 1; // por errores de redondeo
}

void borrarVector(vector<Individuo*>& gen){
    for(vector<int>::size_type j = 0; j != gen.size(); j++) {
        delete gen[j];
    }
    gen.clear();
}

void Generation(vector<Individuo*>& genActual,vector<Individuo*>& genNueva){
    borrarVector(genNueva);
    double sumScaledFitness = 0;
    for(vector<int>::size_type j = 0; j != genActual.size(); j++) {
        sumScaledFitness += genActual[j]->scaled_fitness;
    }
    while(genNueva.size() < MAX_POP){
        int mate1 = Select(genActual,sumScaledFitness);
        int mate2 = Select(genActual,sumScaledFitness);
        Crossover(genActual,mate1,mate2,genNueva);
    }
}

void printDetalle(const vector<Individuo*>& generacion, unsigned int best){
     printf("  idx %20s%10s %7s   %4s  mutac  fitness scaled_fit\n","cromosoma","","dec","padres");
     for(vector<int>::size_type j = 0; j != generacion.size(); j++) {
         Individuo* ind = generacion[j];
         string cromo = CromosomaString(ind);
         if (ind->pos_cross == 0){
             printf(" %c %2d: %s  %.3f     -     %2d  -> %0.5f %0.5f\n",(j == best ? '*' : ' '), (int)j,cromo.c_str(),ind->dec, ind->cant_mutaciones,ind->fitness,ind->scaled_fitness);
         }
         else if (ind->pos_cross == LONG_CROMOSOMA){
             printf(" %c %2d: %s  %.3f  [ %02d ]   %2d  -> %0.5f %0.5f\n",(j == best ? '*' : ' '), (int)j,cromo.c_str(),ind->dec,ind->idx_padre1, ind->cant_mutaciones,ind->fitness,ind->scaled_fitness);
         }
         else{
            printf(" %c %2d: %s  %.3f  %02d-%02d,%02d %2d  -> %0.5f %0.5f\n",(j == best ? '*' : ' '),(int)j,cromo.c_str(),ind->dec,ind->idx_padre1,ind->idx_padre2,ind->pos_cross, ind->cant_mutaciones,ind->fitness,ind->scaled_fitness);
        }
    }
     printf("---------------------------------------------------------------------");
     printf("\n");
}

void getStats(const vector<Individuo*>& generacion, int& best, double& maxFitness,double& avgFitness,double& minFitness){
    double sumFitness = 0;
    maxFitness = 0;
    minFitness = 1000;
    best = 0;
    for(vector<int>::size_type j = 0; j != generacion.size(); j++) {
        double fitness =  generacion[j]->fitness;
        if (fitness > maxFitness){
            maxFitness = fitness;
            best = (int)j;
        }
        if (fitness < minFitness){
            minFitness = fitness;
        }
        sumFitness += fitness;
    }
    avgFitness = sumFitness / generacion.size();

}

void printStats(int nroGen, const vector<Individuo*>& generacion, int best, double maxFitness, double avgFitness, double minFitness){
    string best_cromo = CromosomaString(generacion[best]);
    printf("Gen: %d -> Size: %d Min fit: %f Avg fit: %f MAX fit: %f \n",nroGen,generacion.size(),minFitness,avgFitness,maxFitness);
    printf("   Best: %s (dec = %f)\n",best_cromo.c_str(),generacion[best]->dec);
    printf("\n");

    if (MOSTRAR_DETALLE){
        printDetalle(generacion,best);
    }
}

// dado umax, uavg, umin obtiene valores de a y b
void prescale(double umax, double uavg, double umin, double &a, double &b) {

	double delta;
	if (umin > (fmultiple*uavg - umax) / (fmultiple - 1.0)) {
		delta = umax - uavg;
		a = (fmultiple - 1.0) * uavg / delta;
		b = uavg * (umax - fmultiple*uavg) / delta;
	} else {
		delta = uavg - umin;
		a = uavg / delta;
		b = -umin * uavg / delta;
	}
}


double scale(double u, double a, double b) {
	return (a * u + b);
}

void scalepop(double max, double avg, double min, const vector<Individuo*>& pop) {
	double a;
	double b;
	prescale(max, avg, min, a, b);
	//sumfitness = 0.0;
	for(vector<Individuo*>::size_type j = 0; j < pop.size(); j++) {
		pop[j]->scaled_fitness = scale(pop[j]->fitness, a, b);
		//sumfitness += pop[j]->scaled_fitness;
	}
}

void Correr(double* sumaMaxFit ){
    int pob_inicial = 30;
    int best;
    double avg_fit,min_fit,max_fit;
    //double sumfitness;

    vector<Individuo*> genActual;
    vector<Individuo*> genNueva;
    for(int i = 0; i < pob_inicial; i++){
        Individuo* ind = GenerarRandom();
        genActual.push_back(ind);
    }
    getStats(genActual,best,max_fit,avg_fit,min_fit);
    if (ESCALAR){
        scalepop(max_fit,avg_fit,min_fit,genActual);
    }
    printStats(0,genActual,best,max_fit,avg_fit,min_fit);
    sumaMaxFit[0] += max_fit;
    for(int gen = 1; gen <= MAX_GEN; gen++){
        Generation(genActual,genNueva);

        getStats(genNueva,best,max_fit,avg_fit,min_fit);
        sumaMaxFit[gen] += max_fit;
        if (ESCALAR){
            scalepop(max_fit,avg_fit,min_fit,genNueva);
        }
        printStats(gen,genNueva,best,max_fit,avg_fit,min_fit);

        genActual.swap(genNueva);
    }
}

int main()
{
    /* initialize random seed: */
    srand (time(NULL));
    double maxFits[MAX_GEN + 1] = {}; // inicializa en 0
    int corridas = 200;
    for(int i = 0; i < corridas; i++){
        Correr(maxFits);
    }
    if (ESCALAR){
         printf("Config: Escalar = TRUE, fmultiple = %f\n\n",fmultiple);
    }
    else{
         printf("Config: Escalar = FALSE\n\n");
    }

    printf("Promedio de fitness del mejor de cada gen. sobre %d rounds:\n",corridas );
    for(int i = 0; i <= MAX_GEN; i++){
        printf("   AvgMax gen %2d: %.5f\n",i,maxFits[i] / corridas);
    }
    return 0;

}
