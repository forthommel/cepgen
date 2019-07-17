#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#include <Python.h>

#include "Handler.h"

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  class Limits;
  class ParametersList;
  namespace card
  {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler : public Handler
    {
      public:
        /// Read a standard configuration card
        explicit PythonHandler( const char* file );
        ~PythonHandler();

      private:
        static constexpr const char* MODULE_NAME = "mod_name";

        static constexpr const char* PROCESS_NAME = "process";
        static constexpr const char* HADR_NAME = "hadroniser";
        static constexpr const char* LOGGER_NAME = "logger";
        static constexpr const char* INTEGRATOR_NAME = "integrator";
        static constexpr const char* GENERATOR_NAME = "generator";
        static constexpr const char* OUTPUT_NAME = "output";

        static constexpr const char* PDGLIST_NAME = "PDG";

        static void throwPythonError( const std::string& message );
        static std::string pythonPath( const char* file );
        static PyObject* element( PyObject* obj, const char* key );
        static PyObject* encode( const char* str );

        template<typename T> bool is( PyObject* obj ) const;
        template<typename T> T get( PyObject* obj ) const;
        template<typename T> bool isVector( PyObject* obj ) const;
        template<typename T> std::vector<T> getVector( PyObject* obj ) const;

        void fillParameter( PyObject* parent, const char* key, bool& out );
        void fillParameter( PyObject* parent, const char* key, int& out );
        void fillParameter( PyObject* parent, const char* key, unsigned long& out );
        void fillParameter( PyObject* parent, const char* key, unsigned int& out );
        void fillParameter( PyObject* parent, const char* key, double& out );
        void fillParameter( PyObject* parent, const char* key, std::string& out );
        void fillParameter( PyObject* parent, const char* key, Limits& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<int>& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<double>& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<std::string>& out );
        void fillParameter( PyObject* parent, const char* key, ParametersList& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<ParametersList>& out );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseLogging( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseHadroniser( PyObject* );
        void parseOutputModule( PyObject* );
        void parseExtraParticles( PyObject* );
    };
    template<> bool PythonHandler::is<bool>( PyObject* obj ) const;
    template<> bool PythonHandler::is<int>( PyObject* obj ) const;
    template<> bool PythonHandler::is<long>( PyObject* obj ) const;
    template<> int PythonHandler::get<int>( PyObject* obj ) const;
    template<> unsigned long PythonHandler::get<unsigned long>( PyObject* obj ) const;
    template<> bool PythonHandler::is<ParametersList>( PyObject* obj ) const;
    template<> ParametersList PythonHandler::get<ParametersList>( PyObject* obj ) const;
    template<> bool PythonHandler::is<double>( PyObject* obj ) const;
    template<> double PythonHandler::get<double>( PyObject* obj ) const;
    template<> bool PythonHandler::is<std::string>( PyObject* obj ) const;
    template<> std::string PythonHandler::get<std::string>( PyObject* obj ) const;
    template<> bool PythonHandler::is<Limits>( PyObject* obj ) const;
    template<> Limits PythonHandler::get<Limits>( PyObject* obj ) const;
  }
}

#endif

