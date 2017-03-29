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
const int MAX_GEN = 7;

const double pCross = 0.6;
const double pMutation = 0.01;

class Individuo{
public:
    bool cromosoma[LONG_CROMOSOMA];
    int longitud;
    double fitness;
    double dec;
    int idx_padre1;
    int idx_padre2;
    int pos_cross;
    int cant_mutaciones = 0;

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
int Select(const vector<Individuo*>& pop, double sumfitness ){
    double suma = 0;
    double punto_ruleta = random_float() * sumfitness;

    for(vector<Individuo*>::size_type i = 0; i < pop.size(); i++){
        suma += pop[i]->fitness;
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
    double sumFitness = 0;
    for(vector<int>::size_type j = 0; j != genActual.size(); j++) {
        sumFitness += genActual[j]->fitness;
    }
    while(genNueva.size() < MAX_POP){
        int mate1 = Select(genActual,sumFitness);
        int mate2 = Select(genActual,sumFitness);
        Crossover(genActual,mate1,mate2,genNueva);
    }
}

void printDetalle(const vector<Individuo*>& generacion,int best){
     printf("  idx %20s%10s %7s   %8s    mutac  fitness\n","cromosoma","","dec","padres");
     for(vector<int>::size_type j = 0; j != generacion.size(); j++) {
         Individuo* ind = generacion[j];
         string cromo = CromosomaString(ind);
         if (ind->pos_cross == 0){
             printf(" %c %2d: %s  %f     -     %2d    -> %f\n",(j == best ? '*' : ' '), (int)j,cromo.c_str(),ind->dec, ind->cant_mutaciones,ind->fitness);
         }
         else if (ind->pos_cross == LONG_CROMOSOMA){
             printf(" %c %2d: %s  %f  [ %02d ]   %2d    -> %f\n",(j == best ? '*' : ' '), (int)j,cromo.c_str(),ind->dec,ind->idx_padre1, ind->cant_mutaciones,ind->fitness);
         }
         else{
            printf(" %c %2d: %s  %f  %02d-%02d,%02d %2d    -> %f\n",(j == best ? '*' : ' '),(int)j,cromo.c_str(),ind->dec,ind->idx_padre1,ind->idx_padre2,ind->pos_cross, ind->cant_mutaciones,ind->fitness);
        }
    }
     printf("---------------------------------------------------------------------");
     printf("\n");
}

void printStats(int nroGen, const vector<Individuo*>& generacion){
    double sumFitness = 0;
    double maxFitness = 0;
    int best = 0;
    for(vector<int>::size_type j = 0; j != generacion.size(); j++) {
        double fitness =  generacion[j]->fitness;
        if (fitness > maxFitness){
            maxFitness = fitness;
            best = (int)j;
        }
        sumFitness += fitness;
    }
    double avgFitness = sumFitness / generacion.size();
    string best_cromo = CromosomaString(generacion[best]);
    printf("Gen: %d -> Size: %d Avg fit: %f Max fit: %f\n",nroGen,generacion.size(),avgFitness,maxFitness);
    printf("   Best: %s (%f)\n",best_cromo.c_str(),generacion[best]->dec);
    printf("\n");

    // descomentar para mostrar informacion por cada individuo
    printDetalle(generacion,best);
}

int main()
{
    /* initialize random seed: */
    srand (time(NULL));
    int pob_inicial = 30;
    vector<Individuo*> genActual;
    vector<Individuo*> genNueva;
    for(int i = 0; i < pob_inicial; i++){
        Individuo* ind = GenerarRandom();
        genActual.push_back(ind);
    }
    printStats(0,genActual);

    for(int gen = 1; gen <= MAX_GEN; gen++){
        Generation(genActual,genNueva);
        printStats(gen,genNueva);

        genActual.swap(genNueva);
    }

    return 0;

}
