#ifndef HELPER_H
#define HELPER_H

#include <string_view>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace helper {
  std::vector<std::string> getFileList(std::string_view file_name)
  {
    std::vector<std::string> listOfFiles;

    std::ifstream fileList(file_name.data());

    if (!fileList.is_open()) {
      std::cout << "unable to open file " << file_name << "\n";
      return listOfFiles;
    }

    int numFiles = 0;
    std::string tmp;
    while (std::getline(fileList, tmp)) {
      numFiles++;
    }

    std::cout << "Number of files in the list: " << numFiles << '\n';

    fileList.clear();
    fileList.seekg(0, std::ios::beg);

    listOfFiles.reserve(numFiles);

    std::string line;
    while (std::getline(fileList, line)) {
      listOfFiles.emplace_back(line);
    }

    return listOfFiles;
  }
}

#endif
