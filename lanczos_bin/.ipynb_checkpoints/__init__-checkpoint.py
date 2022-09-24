#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg


from .lanczos import exact_lanczos,lanczos_FA
from .lanczos_RA import streaming_banded_rational, streaming_banded_inv,streaming_banded_prod,streaming_LDL

from .misc import get_discrete_nodes,get_cheb_nodes,model_problem_spectrum,discrete_laplacian_spectrum