/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SCREENOUTPUT_HPP
#define SCREENOUTPUT_HPP

#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

namespace sg {
  namespace base {

    /**
     *  This is used to implement the output on the command line
     *
     *  @version $HEAD$
     */
    class ScreenOutput {
      public:
        /**
         * Standard constructor
         */
        ScreenOutput();

        /**
         * Standard destructor
         */
        ~ScreenOutput();

        /**
         *  update function running windows
         *
         *  @param progress percentage value
         *  @param status status in words
         */
        void update(size_t progress, std::string status);

        /**
         * writes the header of the screen output
         *
         * @param appTitle the application's title
         * @param appAuthor the application's author
         */
        void writeTitle(std::string appTitle, std::string appAuthor);

        /**
         * writes some help information to the console
         *
         * @param helpText the helptext to be displayed
         */
        void writeHelp(std::string helpText);

        /**
         * start the screen output
         *
         * @param text the start text to be displayed
         */
        void writeStartSolve(std::string text);

        /**
         * writes empty lines
         *
         * @param numLines number of empty lines to display
         */
        void writeEmptyLines(size_t numLines);

      private:
#ifdef _WIN32
        /// coordinates of the cusor's position
        COORD pos;
        /// buffers of the CMD
        CONSOLE_SCREEN_BUFFER_INFO info;
        /// handle to the CMD
        HANDLE hCon;
#endif
        /// is this the first run
        bool first_run;
    };

  }
}

#endif /* SCREENOUTPUT_HPP */
