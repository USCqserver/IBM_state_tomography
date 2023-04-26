def get_time(filename):
    vars = filename.split('_')[-7]
    return vars[0:-2]
#convert to integer?

get_time(unsortedQASM[0]) #test if the function works

times = list(set([get_time(x) for x in unsortedQASM])) #remove duplicates
times.sort()
print(times) #check if all the times are covered

import random
random.shuffle(times) #randomize the timesteps
print(times)

sorted_qasm_list = []
for time in times:
    sorted_qasm_list += [folder +  f for f in unsortedQASM if get_time(f) == time]
    #add a sort function for the time clustered groups
    #IC F vs D



tokens = { 
'bbk.pokharel@gmail.com' :  '4c7ba9466e8d1bafc095fbc136f3e59b14e23457ae0883ca7e235347842acfe88ef14342aa83703f8a2def58b54dc0f1b7f963232f88248e6dc547c5e8ca3a65' ,
'pokharel@usc.edu' : '21acf454868c89a86938ffcb6b17cc1236f783c6e340937ab3205598608408623893094770f403dd95ba758d9864ef2135eda8f381804787e15e8ab2aedd3ba3',
# 'Matt 1' : '66a15e358a104a2680592d38599216911e752e373dec894b0f7196d80982e4ed14fd8337d9f10b49f1f4b98bcc7c651b2f2c150726914a3984214097fe70eb94',
# 'Matt 2' : 'a4a855f0dd4e7f3d025718c244e182c9f9aa20d62579ac7b565e371fa36409a5c49ec64e18aa1c54e822e4a714a94a1751998121e068b569b5898c18b7f87f4f',
#     'Matt 3' : '0259e75023ad771da659fe35ccb6efa1c3f179511f31bca63701c4bffe40c2d856ccc092b01bf4d45cf179a10175c45143709cdcf59bdca5b8a3e2e14f3bf6af',
#     'Vinay 1': '31c2289f9d24c098baf2ba21c4534a0811dffaa2c3d350dad246528e0a3f10509253a7047a24e1eaa27fc71ef98954ec7adbbedece7fc422855abcd53a730049',
#     'Vinay 2' : 'be376bd6ecef35c8a887c42528263491dca319feee414e09929b3e60637e529a88d049aef61b835c811f362d21378885afb6b4d3af42c04033b2d6d243c03105',
#     'Vinay 3': '7540e0d1fe55966d0034885ad3b21c25793ae3e66fe71d143803a0db20cc7941fc5d128bba8bbafa5bf59fae83932c0691972aa901141071449d73e50127d59a',
#     'hiyanthi.peiris' : 'e993081e329ba581be46af6c2c739bfdc1c0fbca06081519157f238f8dbb636bfa783adfcbdc141fcd6bc03ed9e675b35d90cc34e8f81b19e3f04066d916e8af',
     'peiris.pokharel' : '56e4e455625f6ad859b235b3821c1cc145124e0b66aaed5d07d07d27856ebe8a76e3f36128dad5c6eae8ad77babac60932cbb8b328c9c95138e70ac69711cc36',
    'haimeng@usc.edu':'bc018d68b58612deebcfb4a7c0597a25f2ecc95a55b00aeb1b70092d643e4a5207079ec12346fa791f2ad8bf58d410305cb82b88531036335d634ef09446fff7',
  #        'zhanghaimeng1994@gmail.com':'96979347f50056d21833fce8259c5b3821c4651248bff09dcd30f4f6364d462fefc479c0389682c214442c637ffc329be381220fd5dd4bf53f4c23bd91d64e39',
#      'theireasychair@gmail.com':'cacd020d84df513852c6f59d4648e7f1e4e4d9c29f032306c60ab23dccf0d225d4c1fa6857c51a0c792e07b2a96d358e6494a23e7db64f20ece90b6eb06f5c85',
 #         'theireasychair@icloud.com':'7a265cd47bf57c072146f15d1967eb165b9704352fc4e4736f2b4f64bb8ff2462be10301cc127739a5a7aa4333e07b52c304e6b452ddc6b15d9c5c2056ef2a13',
 #         'jinweixu@stanford.edu':'6dbc2e530007ad4fcd44abccf232873db3589341cd9d728219c4c97ccc489be727ca3f46412e10a0292f124e895c0e5623471bb1d5c3f54b4eb3bdc9d90e5338'
}
