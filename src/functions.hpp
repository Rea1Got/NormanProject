#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::vector<double> read_input(const std::string filename) {
  std::ifstream inputFile(filename);
  std::vector<double> numbers;
  double num;

  while (inputFile >> num) {
    numbers.push_back(num);
  }
  inputFile.close();

  return numbers;
}
#endif
