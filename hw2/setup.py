from setuptools import extension, setup

setup(
    name= 'mykmeanssp',
    ext_modules=[
        extension(
            'mykmeanssp',
            ['kmeans.c'],
        ),
    ]

)