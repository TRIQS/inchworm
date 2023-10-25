# Generated automatically using the command :
# c++2py ../../c++/inchworm/post_process.hpp --members_read_only -N inchworm -a inchworm -m post_process -o post_process -C triqs --moduledoc="The inchworm postprocess functionality" --cxxflags="-std=c++17" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "post_process", doc = r"The inchworm postprocess functionality", app_name = "inchworm")

# Imports

# Add here all includes
module.add_include("inchworm/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""

using namespace inchworm;
""")




module.generate_code()