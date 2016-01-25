/*
 * JSONException.hpp
 *
 *  Created on: Nov 9, 2015
 *      Author: pfandedd
 */

#pragma once

#include <exception>
#include <string>

#include <sgpp/base/tools/json/Token.hpp>

namespace json {

  class json_exception: public std::exception {
    private:

      std::string message;

    public:
      json_exception(Token& token, const std::string& message);

      json_exception(const std::string& message);

      virtual const char* what() const throw () override;
  };

}
