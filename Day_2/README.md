This is the sample code for Day 2 of the software summer school MM project.

Day 2 covers
- PEP 8
    - make point that all languages have recommended style and formatting.
    - C++ style and formatting - how will this be covered?
- Numpy style docstrings
- Error and Exception handling
- testing using pytest

## Student milestones:
Students will work with their teams on a common repository (MM_teamXX_2019) to fulfill the following milestones. Changes to the repo should be done using a Fork/PR model, where every change must be reviewed by one other person before merging.
1. Write unit tests for your code following the given guidelines:
    - Unit tests which check the values returned from the functions should be written for as many functions as possible. 
        - generate_initial_state
            For each option (`file` and `random`), verify the length of the returned coordinates is what is expected. 
            - For "method=file", test that the first and last line of coordinates are equal to the numbers read by your function.
            - For "method=random", test that the maximum coordinate does not exceed box_length/2 and that the minimum coordinate is not less than -box_length/2. *Hint* - Use the functions numpy.min() and numpy.max() on the coordinates array.
        - lennard_jones_potential
             - Use `@pytest.mark.parametrize` to test the function `lennard_jones_potential` for a range of distances.
                - Hint: Your definition should look something like
                    ~~~
                    @pytest.mark.parametrize("distance2, expected_energy", [
                        ... values
                    ])
                    ~~~
        - minimum_image_distance
            - Use `@pytest.mark.parametrize` to test the function `minimum_image_distance` for a range of distances.
        - total_potential_energy
            - Verify that the total potential energy calculated by your function matches that given by NIST for the sample configuration at a specified cut off.
            - **Challenge** - NIST gives reference energies for multiple cut offs. What pytest structure can you use here to test both cutoffs?
            - **Challenge** - Since both `test_generate_initial_state_file` and `test_potential_energy` read from a reference file, create a @pytest.fixture which reads from this file and returns the coordinates. Use this fixture with the two tests.
        - tail_correction
            - Test that your tail correction matches references from NIST (U_LRC*). https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
            - What pytest feature could you use for this?


1. Add error handling to your functions.
    - For example, in the function `generate_initial_state`, your function should check that input parameters are compatible for each method.
        - if the method is "random", expected inputs are `num_particles` and `box_size`. Having additional or missing arguments should cause a TypeError.
        - if the method is "file", expected inputs are `fname`. Having additional or missing arguments should cause a TypeError.
        - You should check that the method is either "random" or "file". If not either of these, raise a ValueError.
        - use a `try except` clause to open the file for `method=file`
    **Instructor note** - you could walk them through writing this error checking - probably the most complicated of any function.
