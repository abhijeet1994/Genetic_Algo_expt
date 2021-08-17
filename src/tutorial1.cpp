#include "getFoodData.h"
#include "genetic_algo.h"
#include <memory>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problems/schwefel.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/population.hpp>
#include <pagmo/algorithms/nsga2.hpp>

using namespace pagmo;


void interpretRule(const std::vector<std::string> &rule,const FoodDatabase &foodDatabase, std::vector<int> &InterpretedRule);
void debugMaps(const std::map<std::string, int> &foodMap);
std::vector<float> getDecisionVariableConstrain();

struct fitness_database {
    float cal_target;
    float proteins_target = 52.5;
    float carbs_target = 87.5 ;
    float fats_target = 16 ; 
    std::vector<std::vector<foodData>> foodList;
    std::vector<int> ruleIndex; 
    std::vector<std::vector<float>> req_perc = {{0.349, 0.55}, {0.2, 0.349} , {0.2, 0.349}, {0.2, 0.500}};

// Dictionary : Id : Function
// 1  -> Cal_Calc , 2-> carbs ...
// Array  = [1, 2, 3, 7] Refer line 69
// Fitness Function : DBase and Dictionary[ID] from Array call


    vector_double fitness(const vector_double &dv) const
    {
        float cals = 0;
        float proteins = 0;
        float carbs = 0 ;
        float fats = 0 ; 
        float meal_cal ;
        
        int rule_size = ruleIndex.size();
        std::vector<float> cal_ratio ; 
        for (int i = 0 ; i < rule_size ; i ++){
            foodData next_item = foodList[ruleIndex[i]][floorl(dv[2*i])]; 
            float next_serving = (dv[(2 * i) + 1]/100.0);
            meal_cal = (next_item.kiloCalories) * next_serving;
            cals = cals + meal_cal;
            proteins = proteins + ((next_item.proteinContent) * next_serving);
            carbs = carbs + ((next_item.carbsContent) * next_serving);
            fats = fats + ((next_item.FatContent) * next_serving);
            cal_ratio.push_back(meal_cal);
        }

        for (int i = 0 ; i < cal_ratio.size() ; i ++){
            cal_ratio[i] = cal_ratio[i] / cals ; 
            cal_ratio[i] = get_cal_ratio_fitness(cal_ratio[i], req_perc[i]) ; 
        }
        
        
        cals = (cal_target - cals) * (cal_target - cals); 
        proteins = (proteins_target - proteins) * (proteins_target - proteins)  ;
        fats = (fats_target- fats) * (fats_target- fats) ;
        carbs = (carbs_target - carbs) *  (carbs_target - carbs);
        vector_double fit_vec = {cals, proteins, carbs, fats} ; 
        
        //std::cout << fit_vec[0] << " , " << fit_vec[1] << " , " << fit_vec[2] << " , " << fit_vec[3] << " , " << std::endl;
        for (int i = 0 ; i < cal_ratio.size() ; i ++){
            fit_vec.push_back(cal_ratio[i]) ;   
        }
        
        return fit_vec;
    }

    float get_cal_ratio_fitness(const float &a, const std::vector<float> &range) const
    {
        if (a < range[0]){return 5 * (a - range[0]) * (a - range[0]);}
        if (a > range[1]){return 5 * (a - range[1]) * (a - range[1]);}
        float mid = (range[0] + range[1])/2.0;
        return (mid - a) * (mid - a);
    }

};

struct pop_data {
    struct fitness_database fit_db;

    vector_double fitness(const vector_double &dv) const
    {
        return fit_db.fitness(dv);
    }

    void sort_fitness(std::vector<vector_double> &fit_dv)
    {
        vector_double next_ind ; 
        std::vector<std::vector<float>> fit_vec;
        
        for (int i = 0 ; i < fit_dv.size() ; i ++){
            next_ind.clear(); 
            next_ind = fitness(fit_dv[i]);
            std::vector<float> next_vec;
            next_vec.clear() ; 
            for (int j = 0 ; j < next_ind.size() ; j ++){
                next_vec.push_back(float(next_ind[j]));
            }
            fit_vec.push_back(next_vec);
        }

        std::sort(fit_vec.begin(), fit_vec.end(),[&](const std::vector<float> &a, const std::vector<float> &b)
        {
            return a[0] < b[0];
        });


        for (int i = 0 ; i < fit_vec.size() ; i ++) {
            std::vector<float> d = fit_vec[i];
            std::cout << "Ind : " << i ;
            std::cout << ", Cals : " << d[0] + fit_db.cal_target <<", Proteins :" << d[1] + fit_db.proteins_target << ", Carbs : " << d[2] + fit_db.carbs_target <<", Fats: " << d[3] + fit_db.fats_target; 
            
            // << " , Meal_ratio : " ;
            // for (int i = 4 ; i < d.size() ; i ++){
            //     std::cout << d[i] <<" , ";
            // }

            std::cout << "\n";
        }

        // for (int i = 0 ; i < fit_dv.size() ; i ++){
        //     std::vector<float> next_vec;
        //     next_vec.clear();
        //     for (int j = 0 ; j < fit_dv[i].size() ; j ++){
        //         next_vec.push_back(float(fit_dv[i][j]));
        //     }
        //     fit_vec.push_back(next_vec);
        // }
        //for (int i = 0 ; i < fit_vec.size() ; i ++) {
    
    }
};

struct problem_v0 {
    float cal_target;
    std::vector<std::vector<foodData>> foodList;
    std::vector<int> ruleIndex; 
    vector_double fitness(const vector_double &dv) const
    {
        float cals = 0;
        int rule_size = ruleIndex.size();
        for (int i = 0 ; i < rule_size ; i ++){
            foodData next_item = foodList[ruleIndex[i]][floorl(dv[2*i])]; 
            float next_serving = (dv[(2 * i) + 1]/100.0);
            cals = cals + ((next_item.kiloCalories) * next_serving);    
        }
        cals = (cal_target - cals) * (cal_target - cals)  ; 
        if (cals < 0.1){
            cals = 0.1;
        }
        return {cals};    
    }
    
    std::pair<vector_double, vector_double> get_bounds() const
    {
        std::vector<double> lower;
        std::vector<double> upper;
        for (const int &foodCat : ruleIndex){
            lower.push_back(0);
            lower.push_back(0);
            upper.push_back(foodList[foodCat].size() - 1);
            upper.push_back(200);
        }
        return {lower, upper};
    }
};

struct problem_nsga {
    struct fitness_database fit_db;

    vector_double fitness(const vector_double &dv) const
    {
        return fit_db.fitness(dv);
    }

    vector_double::size_type get_nobj() const
    {
        return 8;
    }

    // Implementation of the box bounds.
    std::pair<vector_double, vector_double> get_bounds() const
    {
        std::vector<double> lower;
        std::vector<double> upper;
        for (const int &foodCat : fit_db.ruleIndex){
            lower.push_back(0);
            lower.push_back(0);
            upper.push_back(fit_db.foodList[foodCat].size() - 1);
            upper.push_back(400);
        }
        return {lower, upper};
    }
};

int main()
{
    std::string FileLocation = "/home/abhi/ga-cpp/files/food_data_2.csv";
    getFoodData food_data(FileLocation);
    FoodDatabase foodDB = food_data.getFoodDataList(); 
    std::vector<std::string> rule = {"Main","Veggies","Veggies","Side"};  //.getRule()
    std::vector<int> InterpretedRule;
    InterpretedRule.clear();
    InterpretedRule.reserve(rule.size());
    interpretRule(rule,foodDB,InterpretedRule);
    std::vector<float> decisionVariableConstrain = getDecisionVariableConstrain();

    int POP_SIZE = 1000; 
    float IDEAL_CAL = 700; 
    int GEN = 300;
    
    struct fitness_database fit_data = {} ; 
    fit_data.cal_target = IDEAL_CAL;
    fit_data.proteins_target = 52.5 ;
    fit_data.carbs_target = 87.5 ;
    fit_data.fats_target = 16 ; 
    fit_data.foodList = foodDB.food_database;
    fit_data.ruleIndex = InterpretedRule;

    struct problem_nsga prob = {};
    prob.fit_db = fit_data ;     

    struct pop_data pop_data = {};
    pop_data.fit_db = fit_data ;
    
    problem p{prob};
    algorithm algo{nsga2(GEN, 0.95, 20. , 0.40, 50)};
    algo.set_verbosity(2);
    population pop{p, POP_SIZE, 5u}; 
    pop = algo.evolve(pop);
    
    std::vector<vector_double> comp_pop; 
    for (int i = 0 ; i < POP_SIZE ; i ++) {
        vector_double dv = pop.get_x()[i];
        comp_pop.push_back(dv);
    }
    pop_data.sort_fitness(comp_pop);
}

std::vector<float> getDecisionVariableConstrain()
{
    std::vector<float> decisionVar;
    decisionVar.reserve(4);
    decisionVar[inverseCalDiff] = 1.0F; //
    decisionVar[carbCal] = 1.0F; // 
    decisionVar[proteinCal] = 1.0F; //
    decisionVar[FatCal] = 1.0F; //
    return decisionVar;
}

void interpretRule(const std::vector<std::string> &rule,const FoodDatabase &foodDatabase,std::vector<int> &InterpretedRule)
{
    for (const std::string &ruleComp : rule)
    {
        InterpretedRule.push_back(foodDatabase.foodCategories.at(ruleComp));
    }
}

void debugMaps(const std::map<std::string, int> &foodMap)
{
    for (auto const&x : foodMap)
    {
        std::cout << x.first << ":" << x.second <<std::endl;
    }
}