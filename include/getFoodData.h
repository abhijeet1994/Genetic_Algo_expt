#ifndef GETFOODDATA_H
#define GETFOODDATA_H
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <map>

namespace foodDataConsts
{
    const uint NumOfCategories = 20;
    const uint varitiesOfFood = 50;
    const uint NumOfColumInFoodData = 15;
}

struct foodData
{
    int serialNumber; //0
    std::string ingredientName; //1
    std::string classification; //Main,side,Veggies etc //2
    std::string type; //Veg, NonVeg, Vegan //3
    float proteinContent; //4
    float carbsContent; //5
    float FatContent; //6
    float kiloCalories; //7
    float price; //8
    float quantity; //9
};

struct FoodDatabase
{
    std::map<std::string, int> foodCategories;
    std::vector<std::vector<foodData>> food_database;
};
class getFoodData
{
    public:
    getFoodData(std::string &aFileLocation);
    FoodDatabase getFoodDataList();


private:
    FoodDatabase _mFoodData;
    std::string _mFileLocation;
    std::string _mDlimeter =",";
    void getCSVdata();
    foodData storeDataInVector(std::vector<std::string> &inputvector);
};

#endif