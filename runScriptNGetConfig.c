#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "runScriptNGetConfig.h"

config_parameters_t *get_config_parameters(void)
{
    config_parameters_t *cParms;
    cParms = (config_parameters_t *)malloc(sizeof(config_parameters_t));
    
    cParms->filter_respLen = 25;
    cParms->filter_SavitzkyGolay_poly_order = 6;
    cParms->filter_SavitzkyGolay_derivative_degree = 1;
    
    return cParms;
}

int free_config_parameters(config_parameters_t *cParms)
{
    free(cParms);
    return 0;
}

