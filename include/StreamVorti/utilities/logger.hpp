/*
 * StreamVorti - Software for solving PDEs using explicit methods.
 * Copyright (C) 2017  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors (alphabetically):
 *      George C. BOURANTAS
 *      Konstantinos A. MOUNTRIS
 */

/*!
   \file logger.hpp
   \brief Logger class header file.
   \author Konstantinos A. Mountris
   \date 12/01/2018
*/

#ifndef STREAMVORTI_UTILITIES_LOGGER_HPP_
#define STREAMVORTI_UTILITIES_LOGGER_HPP_

#include <string>

namespace StreamVorti {

/*!
 *  \addtogroup Utilities
 *  @{
 */


/*!
 * \class Logger
 * \brief Class implemmenting output messages for logging of StreamVorti.
 */

class Logger
{
public:

    /*!
     * \brief Logger default constructor.
     */
    Logger() {}


    /*!
     * \brief Logger default destructor.
     */
    virtual ~Logger() {}


    /*!
     * \brief Log an StreamVorti message.
     * \param [in] msg The message to be logged.
     * \return [std::string] The logged message.
     */
    inline static std::string Message(const std::string &msg) {
        return "[StreamVorti] " + msg;
    }


    /*!
     * \brief Log an StreamVorti error.
     * \param [in] err The error to be logged.
     * \return [std::string] The logged error.
     */
    inline static std::string Error(const std::string &err) {
        return "[StreamVorti ERROR] " + err;
    }


    /*!
     * \brief Log an StreamVorti warning.
     * \param [in] wrng The warning to be logged.
     * \return [std::string] The logged warning.
     */
    inline static std::string Warning(const std::string &wrng) {
        return "[StreamVorti WARNING] " + wrng;
    }


};

/*! @} End of Doxygen Groups*/

} //end of namespace StreamVorti

#endif //STREAMVORTI_UTILITIES_LOGGER_HPP_
